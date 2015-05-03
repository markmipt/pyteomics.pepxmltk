#!/bin/sh

ver=$(cat VERSION)
echo "Generating source distribution ..."
python setup.py sdist
tarname="pyteomics.pepxmltk-${ver}.tar.gz"
ln -s "dist/${tarname}" "$tarname"
echo "Calculating MD5 sum ..."
md5=$(md5sum "$tarname" | cut -d' ' -f1)
sed -i "s/^md5sums=.*/md5sums=('${md5}')/" PKGBUILD
echo "Generating Python 3 AUR ball ..."
sed -i "s/^pkgver=.*/pkgver=${ver}/" PKGBUILD
mkaurball -f
echo "Patching PKGBUILD ..."
sed -i.old -E "/^source/!s/python(\b)/python2\1/g" PKGBUILD
sed -i.old -E "/^source/!s/python(\b)/python2\1/g" pepxmltk.install
echo "Generating Python 2 AUR ball ..."
mkaurball -f
echo "Restoring PKGBUILD ..."
mv PKGBUILD.old PKGBUILD
mv pepxmltk.install.old pepxmltk.install
rm "$tarname"
