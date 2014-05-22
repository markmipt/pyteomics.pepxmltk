import unittest
from pyteomics.tandem2xml import *
from pyteomics import parser


class ProteaseTest(unittest.TestCase):
    def setUp(self):
        self.protease1 = Protease('[RK]|{P}')
        self.protease2 = Protease('{P}|[KR]')
        self.protease3 = Protease('[RK]|[X]')
        self.protease4 = Protease('[RK]|{X}')
        self.protease5 = Protease('[X]|[X]')

    def test_simple_cut(self):
        self.assertEqual((set('RK'), set('P'), 'C'),
                         (set(self.protease1.cut), set(self.protease1.no_cut), self.protease1.sense))

    def test_n_term_cut(self):
        self.assertEqual((set('RK'), set('P'), 'N'),
                         (set(self.protease2.cut), set(self.protease2.no_cut), self.protease2.sense))

    def test_x_letters(self):
        self.assertEqual((set('RK'), set(''), 'C'),
                         (set(self.protease3.cut), set(self.protease3.no_cut), self.protease3.sense))
        self.assertEqual((set('RK'), set(parser.std_amino_acids), 'C'),
                         (set(self.protease4.cut), set(self.protease4.no_cut), self.protease4.sense))
    def test_cleavage_all(self):
        self.assertEqual((set(parser.std_amino_acids), set(''), 'C'),
                         (set(self.protease5.cut), set(self.protease5.no_cut), self.protease5.sense))

if __name__ == '__main__':
    unittest.main()

