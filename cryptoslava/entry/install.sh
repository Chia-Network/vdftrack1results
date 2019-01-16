#!/bin/sh
set -e

sudo apt-get install -y libmpfr-dev autoconf m4 libtool-bin binutils bison
VDF_BASE="$(pwd)"
VDF_PREFIX="$(pwd)/vdfroot"
if [ -z "$ARCH_CFLAGS" ]; then
	ARCH_CFLAGS='-march=skylake -mtune=skylake'
fi
export CFLAGS="-O3 -g $ARCH_CFLAGS"

(cd gmp && (autoreconf -i || autoreconf -i) && ./configure --disable-shared --prefix=$VDF_PREFIX && PATH=$VDF_BASE:$PATH make -j4 && PATH=$VDF_BASE:$PATH make install)
(cd flint2 && ./configure --disable-shared --prefix=$VDF_PREFIX --with-gmp=$VDF_PREFIX && make -j4 && make install)

export CFLAGS="$CFLAGS -I$VDF_PREFIX/include"
export LDFLAGS="-L$VDF_PREFIX/lib"
make -C entry

echo "--- END OF install.sh ---"
echo "--- END OF install.sh (stderr) ---" >&2
