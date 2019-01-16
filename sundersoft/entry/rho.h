void test_form(integer d) {
    form f=form::generator(d);
    for (int x=0;x<4;++x) {
        f=square(f); //g^2, g^4, g^8, g^16
    }

    form f2=form::identity(d);
    for (int x=0;x<16;++x) {
        f2=multiply(f2, form::generator(d)); //g, g^2, ..., g^16
    }

    assert(f==f2);

    form f3=f.inverse();
    form f4=f3*f;
    form f5=f*f3;
    assert(f4==form::identity(d));
    assert(f5==form::identity(d));
}

struct generator_table {
    integer d;
    vector<form> values;

    generator_table(const integer& t_d, int num_bits) {
        d=t_d;
        form c=form::generator(d);
        for (int x=0;x<num_bits;++x) {
            values.push_back(c);
            c=c*c;
        }
    }
};

struct exponent {
    integer bits;
    form value;
    const generator_table* c_table=nullptr;

    exponent() {}

    exponent(const generator_table& t_table) : c_table(&t_table) {
        bits=0;
        value=form::identity(c_table->d);
    }

    bool get_bit(int index) {
        return bits.get_bit(index);
    }

    void set_bit(int index, bool new_value) {
        bool old_value=get_bit(index);

        if (old_value==new_value) {
            return;
        }

        form f=c_table->values.at(index);
        if (!new_value) {
            f=f.inverse();
        }
        value=value*f;

        bits.set_bit(index, new_value);
    }

    void add_bit(int index) {
        const form& f=c_table->values.at(index);

        //form new_value;
        form new_value=value*f;
        bits+=integer(1)<<index;

        if (test_gpu) {
            gpu_form<half_d_num_bits> value_gpu=value;
            gpu_form<half_d_num_bits> f_gpu=f;
            gpu_form<half_d_num_bits> new_value_gpu=multiply(value_gpu, f_gpu);
            assert(form(new_value_gpu)==new_value);
            //new_value=new_value_gpu;
        }

        value=new_value;
    }

    void set_bits(integer t_bits) {
        *this=exponent(*c_table);
        for (int x=0;x<t_bits.num_bits();++x) {
            if (t_bits.get_bit(x)) {
                set_bit(x, true);
            }
        }
    }

    void assert_valid() {
        exponent expected(*c_table);
        expected.set_bits(bits);

        assert(expected.value==value);
        value.assert_valid(c_table->d);
    }
};

struct thread_state {
    exponent state;
    int thread_number=-1;
    int64 sequence_number=0;

    bool operator<(const thread_state& s) const {
        return state.value<s.state.value;
    }

    bool operator==(const thread_state& s) const {
        return state.value==s.state.value;
    }

    bool is_distinguished(int distinguished_num_bits) const {
        int new_hash=state.value.hash();
        return (new_hash & ((1<<distinguished_num_bits)-1)) == 0;
    }
};

//this is the only thing that runs on the gpu
struct rho_thread {
    thread_state state;
    vector<thread_state> distinguished_states;

    int exponent_num_bits=-1;
    int distinguished_num_bits=-1;

    rho_thread(generator_table& t_table, int t_exponent_num_bits, int t_distinguished_num_bits, int t_thread_number) {
        exponent_num_bits=t_exponent_num_bits;
        distinguished_num_bits=t_distinguished_num_bits;
        state.thread_number=t_thread_number;
        state.state=exponent(t_table);
    }

    void init() {
        state.state.set_bits(rand_integer(exponent_num_bits));
        if (state.is_distinguished(distinguished_num_bits)) {
            distinguished_states.push_back(state);
        }
    }

    void advance() {
        int c_hash=state.state.value.hash();
        int add_index=c_hash%exponent_num_bits;
        state.state.add_bit(add_index);
        ++state.sequence_number;

        if (state.is_distinguished(distinguished_num_bits)) {
            distinguished_states.push_back(state);
        }
    }
};

struct database {
    generator_table c_table;
    int exponent_num_bits=-1;
    int distinguished_num_bits=-1;
    int64 checkpoint_interval=-1;

    //int64 total_num_checkpoints=0;
    int64 checkpoint_sequence_number=0;
    vector<thread_state> checkpoint;

    multiset<thread_state> distinguished_states;
    vector<vector<multiset<thread_state>::iterator>> distinguished_states_by_thread;
    vector<array<thread_state, 2>> collisions;

    //these are in the same state as another thread but are behind it
    set<int> broken_threads;

    database(
        integer t_discriminant,
        int t_exponent_num_bits, int t_table_num_bits, int t_distinguished_num_bits, int64 t_checkpoint_interval
    ) :
        c_table(t_discriminant, t_table_num_bits)
    {
        exponent_num_bits=t_exponent_num_bits;
        distinguished_num_bits=t_distinguished_num_bits;
        checkpoint_interval=t_checkpoint_interval;
        checkpoint_sequence_number=-checkpoint_interval;
    }

    rho_thread make_thread(int thread_number) {
        return rho_thread(c_table, exponent_num_bits, distinguished_num_bits, thread_number);
    }

    void add_distinguished_state(thread_state& s) {
        s.state.assert_valid();
        assert(s.is_distinguished(distinguished_num_bits));

        auto i=distinguished_states.insert(s);

        while (distinguished_states_by_thread.size()<=s.thread_number) {
            distinguished_states_by_thread.emplace_back();
        }

        {
            auto& v=distinguished_states_by_thread.at(s.thread_number);
            assert(v.empty() || v.back()->sequence_number<s.sequence_number);
            v.push_back(i);
        }

        auto i2=i;
        if (i2!=distinguished_states.begin()) {
            --i2;
            if (*i2==s) {
                collisions.push_back({*i2, s});
            }
        }

        auto i3=i;
        ++i3;
        if (i3!=distinguished_states.end()) {
            if (*i3==s) {
                collisions.push_back({*i3, s});
            }
        }
    }

    void add_checkpoint(vector<rho_thread>& threads) {
        vector<thread_state> c;
        for (int x=0;x<threads.size();++x) {
            rho_thread& t=threads[x];

            assert(t.exponent_num_bits==exponent_num_bits);
            assert(t.distinguished_num_bits==distinguished_num_bits);
            //assert(t.state.sequence_number==total_num_checkpoints*checkpoint_interval);
            assert(t.state.sequence_number==checkpoint_sequence_number+checkpoint_interval);
            assert(t.state.thread_number==x);
            t.state.state.assert_valid(); //do this in the cpu to verify that the gpu code works

            for (thread_state& s : t.distinguished_states) {
                assert(s.sequence_number>checkpoint_sequence_number);
                assert(s.sequence_number<=checkpoint_sequence_number+checkpoint_interval);
                //assert(s.sequence_number>(total_num_checkpoints-1)*checkpoint_interval);
                //assert(s.sequence_number<=total_num_checkpoints*checkpoint_interval);

                //broken_threads isn't persisted
                //if the program is restarted from the last checkpoint while there are broken threads, the broken threads will be replayed
                // again before being detected as broken and added to the broken threads set
                if (!broken_threads.count(s.thread_number)) {
                    add_distinguished_state(s);
                }
            }
            t.distinguished_states.clear();

            c.push_back(t.state);
        }

        checkpoint=c;
        //++total_num_checkpoints;

        checkpoint_sequence_number+=checkpoint_interval;
    }

    bool process_potential_collision(const thread_state& a, const thread_state& b) {
        assert(a.state.value==b.state.value);

        if (a.state.bits==b.state.bits) {
            return false;
        }

        // g^a = g^b
        // g^a*g^(-b) = 1
        // g^(a-b) = 1
        //
        // g^a = g^b
        // 1 = g^b*g^(-a)
        // 1 = g^(b-a)
        //
        // g^(|a-b|) = 1

        integer final_bits=(a.state.bits>b.state.bits)? a.state.bits-b.state.bits : b.state.bits-a.state.bits;
        assert(final_bits!=integer(0));

        exponent final_check(c_table);
        final_check.set_bits(final_bits);
        final_check.assert_valid();
        assert(final_check.value==form::identity(c_table.d));

        cout << c_table.values.at(0).a.to_string_dec() << " ";
        cout << c_table.values.at(0).b.to_string_dec() << " ";
        cout << c_table.values.at(0).c.to_string_dec() << " ";
        cout << final_bits.to_string_dec() << "\n";

        return true;
    }

    bool process_collisions_impl(thread_state& s, multiset<thread_state>& all_states) {
        //s is a distinguished state
        //roll s back to the previous distinguished state, and iterate all of the states between it and s's current state
        //add them to all_states and check for a collision. there should be one

        auto& states=distinguished_states_by_thread.at(s.thread_number);
        int s_pos=-1;

        for (int x=0;x<states.size();++x) {
            assert(states[x]->thread_number==s.thread_number);
            if (states[x]->sequence_number==s.sequence_number) {
                assert(*states[x]==s);
                s_pos=x;
                break;
            }
        }

        assert(s_pos!=-1);
        if (s_pos==0) {
            //can't go backwards if it is 0
            print( "Collision ignored because s_pos==0." );
            return false;
        }

        rho_thread s_thread=make_thread(s.thread_number);
        s_thread.state=*states.at(s_pos-1); //distinguished state before s

        assert(s_thread.state.sequence_number<s.sequence_number);

        int64 num_iterations=(s.sequence_number-s_thread.state.sequence_number);
        if (max_process_collision_iterations!=0 && num_iterations>max_process_collision_iterations) {
            print(
                "Collision ignored because there are too many iterations between distinguished points. num_iterations=", num_iterations
            );
            return false;
        }

        print(
            "Replaying:",
            "Thread number: ", s_thread.state.thread_number,
            "Start sequence number: ", s_thread.state.sequence_number,
            "End sequence number: ", s.sequence_number
        );

        while (s_thread.state.sequence_number<s.sequence_number) {
            //can't have a collision at the first state or it would have gotten detected already
            s_thread.advance();
            auto i=all_states.find(s_thread.state);
            if (i!=all_states.end() && process_potential_collision(s_thread.state, *i)) {
                return true;
            }
            all_states.insert(s_thread.state);
        }
        assert(s_thread.state==s);

        return false;
    }

    bool process_collisions() {
        while (!collisions.empty()) {
            array<thread_state, 2> collision;
            int64 collision_sequence_number=-1;
            int collision_pos=-1;

            //want the earliest collision if there are multiple collisions
            for (int x=0;x<collisions.size();++x) {
                auto& c=collisions[x];
                int64 c_sequence_number=min(c[0].sequence_number, c[1].sequence_number);
                if (collision_sequence_number==-1 || collision_sequence_number>c_sequence_number) {
                    collision=c;
                    collision_sequence_number=c_sequence_number;
                    collision_pos=x;
                }
            }
            assert(collision_pos!=-1);
            collisions.erase(collisions.begin()+collision_pos);

            assert(collision[0].state.value==collision[1].state.value); //the bits could be different but unlikely

            bool is_broken=
                broken_threads.count(collision[0].thread_number) ||
                broken_threads.count(collision[1].thread_number)
            ;
            if (is_broken) {
                continue;
            }

            print( "Processing collision..." );

            multiset<thread_state> all_states;
            if (process_collisions_impl(collision[0], all_states)) {
                return true;
            }

            if (process_collisions_impl(collision[1], all_states)) {
                return true; //likely to happen
            }

            //if the collision wasn't able to be processed, one of the threads is now useless. this is fine for the gpu since it has
            // lots of threads so it isn't handled
            //could reinitialize the thread to a new random value if i wanted to

            //will disable whatever thread is farther behind to avoid having to replay it all of the time
            //if a thread got into a loop and the collision couldn't be used, it is also broken
            if (collision[0].sequence_number<collision[1].sequence_number) {
                broken_threads.insert(collision[0].thread_number);
            } else {
                broken_threads.insert(collision[1].thread_number);
            }

            print(
                "Useless collision",
                "Thread 0:", collision[0].thread_number,
                "Sequence number 0:", collision[0].sequence_number,
                "Thread 1:", collision[1].thread_number,
                "Sequence number 1:", collision[1].sequence_number,
                "Num broken threads:", broken_threads.size()
            );

            //can't make any progress if all of the threads are broken
            assert(broken_threads.size()<num_threads);
        }

        return false;
    }

    void save_thread_state(ostream& out, const thread_state& c) {
        out << c.state.bits.to_string_dec() << "\n";
        out << c.state.value.a.to_string_dec() << "\n";
        out << c.state.value.b.to_string_dec() << "\n";
        out << c.state.value.c.to_string_dec() << "\n";
        out << c.thread_number << "\n";
        out << c.sequence_number << "\n";
    }

    template<class type> type get_stream(istream& in) {
        type res;
        in >> res;
        return res;
    }

    thread_state load_thread_state(istream& in) {
        thread_state res;
        res.state.bits=integer(get_stream<string>(in));
        res.state.value.a=integer(get_stream<string>(in));
        res.state.value.b=integer(get_stream<string>(in));
        res.state.value.c=integer(get_stream<string>(in));
        res.state.c_table=&c_table;
        res.thread_number=get_stream<int>(in);
        res.sequence_number=get_stream<int>(in);

        //res.state.assert_valid(); //slow; already asserted this anyway
        return res;
    }

    void save() {
        string file_name =
            "rho_database_" + (-c_table.d).to_string_dec() + "_" +
            to_string(checkpoint_sequence_number + 1000000000000ll) + ".log"
        ;
        ofstream out(file_name);
        print( "Saving database to", file_name );

        //out << total_num_checkpoints << "\n";
        out << 0 << "\n";
        out << checkpoint.size() << "\n";
        for (thread_state& c : checkpoint) {
            save_thread_state(out, c);
        }

        out << distinguished_states_by_thread.size() << "\n";
        for (auto& c : distinguished_states_by_thread) {
            out << c.size() << "\n";
            for (auto& d : c) {
                save_thread_state(out, *d);
            }
        }

        assert(collisions.empty());
    }

    void load(ifstream& in) {
        assert(checkpoint.empty());
        assert(collisions.empty());
        assert(distinguished_states_by_thread.empty());
        assert(distinguished_states.empty());

        //total_num_checkpoints=get_stream<int>(in);
        get_stream<int>(in);

        checkpoint_sequence_number=-1;

        int checkpoint_size=get_stream<int>(in);
        for (int x=0;x<checkpoint_size;++x) {
            checkpoint.push_back(load_thread_state(in));

            int64 sequence_number=checkpoint.back().sequence_number;
            assert(sequence_number==checkpoint_sequence_number || checkpoint_sequence_number==-1);
            checkpoint_sequence_number=sequence_number;
        }

        int distinguished_states_by_thread_size=get_stream<int>(in);
        for (int x=0;x<distinguished_states_by_thread_size;++x) {
            distinguished_states_by_thread.emplace_back();

            int c_size=get_stream<int>(in);
            for (int y=0;y<c_size;++y) {
                auto i=distinguished_states.insert(load_thread_state(in));
                distinguished_states_by_thread.back().push_back(i);
            }
        }
    }

    vector<rho_thread> restore_checkpoint() {
        vector<rho_thread> res;
        for (thread_state& c : checkpoint) {
            rho_thread t=make_thread(c.thread_number);
            t.state=c;
            res.push_back(t);
        }

        return res;
    }
};