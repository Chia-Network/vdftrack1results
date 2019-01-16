// Copyright 2018 POA Networks Ltd.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#[macro_use]
extern crate criterion;

use classgroup::{gmp_classgroup::GmpClassGroup, ClassGroup};
use criterion::Criterion;
use gmp::mpz::Mpz;
use std::str::FromStr;

fn bench_square(c: &mut Criterion) {
    for _ in 0..2 {
        let m_2048 = -Mpz::from_str(
            "201493927071865251625903550712920535753645598483515670853547009\
             878440933309489362800393797428711071833308081461824159206915864\
             150805748296170245037221957772328044276705571745811271212292422\
             075849739248257870371300001313586036515879618764093772248760562\
             386804073478433157526816295216137723803793411828867470089409596\
             238958950007370719325959579892866588928887249912429688364409867\
             895510817680171869190054122881274299350947669820596157115994418\
             034091728887584373727555384075665624624856766441009974642693066\
             751400054217209981490667208950669417773785631693879782993019167\
             69407006303085854796535778826115224633447713584423",
        )
        .unwrap();

        let m_1024 = -Mpz::from_str(
            "-11208471744389096429663063172516742066731683613191418514476174383781\
             682509882427394963852743081347678693241523614532942268295868231081182\
             819214054220080323345750407342623884342617809879459211722505867733607\
             400509994975706778681543998242335468203860240586171413971485860382901\
             6409314686266660248501773529803183",
        )
        .unwrap();
        let group_1024 = GmpClassGroup::generator_for_discriminant(m_1024);
        let group_2048 = GmpClassGroup::generator_for_discriminant(m_2048);
        let (group_1024_clone, group_2048_clone) = (group_1024.clone(), group_2048.clone());
        c.bench_function("square 1024", move |b| {
            b.iter(|| group_1024_clone.clone().square())
        });
        c.bench_function("multiply 1024", move |b| {
            b.iter(|| &group_1024 * &group_1024)
        });
        c.bench_function("square 2048", move |b| {
            b.iter(|| group_2048_clone.clone().square())
        });
        c.bench_function("multiply 2048", move |b| {
            b.iter(|| &group_2048 * &group_2048)
        });
    }
}

criterion_group!(benches, bench_square);
criterion_main!(benches);
