# Maintainer: Lev Levitsky <levlev at mail dot ru>
pkgname="python-pyteomics.pepxmltk"
pkgver=0.2.1
pkgrel=1
pkgdesc='Convert X!Tandem XML files to pepXML, create pepXML from Python, run X!Tandem from Python and command line'
arch=('any')
url="https://pypi.python.org/pypi/pyteomics.pepxmltk"
license=('Apache')
depends=('python-pyteomics')
optdepends=('xtandem: the X!Tandem search engine')
options=(!emptydirs)
source=("https://pypi.python.org/packages/source/p/pyteomics/pyteomics.pepxmltk-${pkgver}.tar.gz")
md5sums=('c223a3d02f4762e112634b49e1a80601')
package() {
  cd "${srcdir}/pyteomics.pepxmltk-${pkgver}"
  python setup.py install --root="$pkgdir/" --optimize=1
}

# vim:set ts=2 sw=2 et:
