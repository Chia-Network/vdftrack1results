#!/bin/bash

./gp-sta -q <<HEREDOC
{
d=$1;
n=$2;
\\\\d=-0xdc2a335cd2b355c99d3d8d92850122b3d8fe20d0f5360e7aaaecb448960d57bcddfee12a229bbd8d370feda5a17466fc725158ebb78a2a7d37d0a226d89b54434db9c3be9a9bb6ba2c2cd079221d873a17933ceb81a37b0665b9b7e247e8df66bdd45eb15ada12326db01e26c861adf0233666c01dec92bbb547df7369aed3b1fbdff867cfc670511cc270964fbd98e5c55fbe0947ac2b9803acbfd935f3abb8d9be6f938aa4b4cc6203f53c928a979a2f18a1ff501b2587a93e95a428a107545e451f0ac6c7f520a7e99bf77336b1659a2cb3dd1b60e0c6fcfffc05f74cfa763a1d0af7de9994b6e35a9682c4543ae991b3a39839230ef84dae63e88d90f457;
\\\\n=2097152;

p=-d;
L=sqrtint(sqrtint(p+1)+1)+1;
x=Qfb(2,1,(p+1)/8);

for(i=1,n,
  x = qfbnucomp(x,x,L);
);

v = Vec(x);
print(v[1]," ",v[2]);

quit;
}

HEREDOC

