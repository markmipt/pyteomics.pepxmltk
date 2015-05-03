# Maintainer: Lev Levitsky <levlev at mail dot ru>
pkgname="python-pyteomics.pepxmltk"
pkgver=0.2.2
pkgrel=1
pkgdesc='Convert X!Tandem XML files to pepXML, create pepXML from Python, run X!Tandem from Python and command line'
arch=('any')
url="https://pypi.python.org/pypi/pyteomics.pepxmltk"
license=('Apache')
depends=('python')
optdepends=('xtandem: the X!Tandem search engine (for runtandem)'
            'python-pyteomics: for pepxmltk functionality'
            'python-jinja: for pepxmltk functionality')
options=(!emptydirs)
source=("https://pypi.python.org/packages/source/p/pyteomics.pepxmltk/pyteomics.pepxmltk-${pkgver}.tar.gz")
install=pepxmltk.install
md5sums=('cb691d326be5b0163051580bf14e1d1e')
package() {
  cd "${srcdir}/pyteomics.pepxmltk-${pkgver}"
  python setup.py install --root="$pkgdir/" --optimize=1
}

# vim:set ts=2 sw=2 et:
