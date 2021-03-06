import Ficha4
import unittest

class TestFicha4(unittest.TestCase):

    def test_ler_seq(self):
        self.assertEqual(Ficha4.ler_seq('sequence.txt'), 'ATGGCAGATTCCAACGGTACTATTACCGTTGAAGAGCTTAAAAAGCTCCTTGAACAATGGAACCTAGTAATAGGTTTCCTATTCCTTACATGGATTTGTCTTCTACAATTTGCCTATGCCAACAGGAATAGGTTTTTGTATATAATTAAGTTAATTTTCCTCTGGCTGTTATGGCCAGTAACTTTAGCTTGTTTTGTGCTTGCTGCTGTTTACAGAATAAATTGGATCACCGGTGGAATTGCTATCGCAATGGCTTGTCTTGTAGGCTTGATGTGGCTCAGCTACTTCATTGCTTCTTTCAGACTGTTTGCGCGTACGCGTTCCATGTGGTCATTCAATCCAGAAACTAACATTCTTCTCAACGTGCCACTCCATGGCACTATTCTGACCAGACCGCTTCTAGAAAGTGAACTCGTAATCGGAGCTGTGATCCTTCGTGGACATCTTCGTATTGCTGGACACCATCTAGGACGCTGTGACATCAAGGACCTGCCTAAAGAAATCACTGTTGCTACATCACGAACGCTTTCTTATTACAAATTGGGAGCTTCGCAGCGTGTAGCAGGTGACTCAGGTTTTGCTGCATACAGTCGCTACAGGATTGGCAACTATAAATTAAACACAGACCATTCCAGTAGCAGTGACAATATTGCTTTGCTTGTACAGTAA')
        self.assertEqual(Ficha4.ler_seq('.txt'), 'O ficheiro não pode ser aberto')
        self.assertEqual(Ficha4.ler_seq(''), 'Introduza um ficheiro válido')

    def test_ler_FASTA_seq(self):
        self.assertEqual(Ficha4.ler_FASTA_seq('sequence.fasta'), 'ATGGCAGATTCCAACGGTACTATTACCGTTGAAGAGCTTAAAAAGCTCCTTGAACAATGGAACCTAGTAATAGGTTTCCTATTCCTTACATGGATTTGTCTTCTACAATTTGCCTATGCCAACAGGAATAGGTTTTTGTATATAATTAAGTTAATTTTCCTCTGGCTGTTATGGCCAGTAACTTTAGCTTGTTTTGTGCTTGCTGCTGTTTACAGAATAAATTGGATCACCGGTGGAATTGCTATCGCAATGGCTTGTCTTGTAGGCTTGATGTGGCTCAGCTACTTCATTGCTTCTTTCAGACTGTTTGCGCGTACGCGTTCCATGTGGTCATTCAATCCAGAAACTAACATTCTTCTCAACGTGCCACTCCATGGCACTATTCTGACCAGACCGCTTCTAGAAAGTGAACTCGTAATCGGAGCTGTGATCCTTCGTGGACATCTTCGTATTGCTGGACACCATCTAGGACGCTGTGACATCAAGGACCTGCCTAAAGAAATCACTGTTGCTACATCACGAACGCTTTCTTATTACAAATTGGGAGCTTCGCAGCGTGTAGCAGGTGACTCAGGTTTTGCTGCATACAGTCGCTACAGGATTGGCAACTATAAATTAAACACAGACCATTCCAGTAGCAGTGACAATATTGCTTTGCTTGTACAGTAA')
        self.assertEqual(Ficha4.ler_FASTA_seq('.txt'), 'Este ficheiro não corresponde a um formato FASTA')
        self.assertEqual(Ficha4.ler_FASTA_seq(''), 'Introduza um ficheiro no formato FASTA')

    def complemento_inverso(self):
        self.assertEqual(Ficha4.complemento_inverso('ATGCATAGGTCGACCGATTGA'), 'TCAATCGGTCGACCTATGCAT')
        self.assertEqual(Ficha4.complemento_inverso('atgcataggtcgaccgattga'), 'TCAATCGGTCGACCTATGCAT')
        self.assertEqual(Ficha4.complemento_inverso('heufbe23023r231fj1n1'), 'Insira uma sequência de DNA válida')

    def transcricao(self):
        self.assertEqual(Ficha4.transcricao('ATGCATAGGTCGACCGATTGA'), 'AUGCAUAGGUCGACCGAUUGA')
        self.assertEqual(Ficha4.transcricao('atgcataggtcgaccgattga'), 'AUGCAUAGGUCGACCGAUUGA')
        self.assertEqual(Ficha4.transcricao('3214FHDF9Q63HE9sdshd'), 'Insira uma sequência de DNA válida')

    def traducao(self):
        self.assertEqual(Ficha4.traducao('ATGCATAGGTCGACCGATTGA'), 'MHRSTD_')
        self.assertEqual(Ficha4.traducao('atgcataggtcgaccgattga'), 'MHRSTD_')
        self.assertEqual(Ficha4.traducao('wjd712312ijefd8f71md'), 'Introduza uma sequência de DNA válida')

    def validar(self):
        self.assertEqual(Ficha4.validar('ATGCATAGGTCGACCGATTGA'), True)
        self.assertEqual(Ficha4.validar('ohduiyfuvdimocpqweiuq'), False)
        self.assertEqual(Ficha4.validar('319123yscnw173t12cnsc'), False)
        self.assertEqual(Ficha4.validar(''), False)
        
    def contar_bases(self):
        self.assertEqual(Ficha4.contar_bases('ATGCATAGGTCGACCGATTGA'), {'A': 6, 'T': 5, 'G': 6, 'C': 4})
        self.assertEqual(Ficha4.contar_bases('atgcataggtcgaccgattga'), {'A': 6, 'T': 5, 'G': 6, 'C': 4})
        self.assertEqual(Ficha4.contar_bases('AUGCAUAGGUCGACCGAUUGA'), {'A': 6, 'U': 5, 'G': 6, 'C': 4})
        self.assertEqual(Ficha4.contar_bases('augcauaggucgaccgauuga'), {'A': 6, 'U': 5, 'G': 6, 'C': 4})
        self.assertEqual(Ficha4.contar_bases('MHRSTD'), {'M': 1, 'H': 1, 'R': 1, 'S': 1, 'T': 1, 'D': 1})
        self.assertEqual(Ficha4.contar_bases('mhrstd'), {'M': 1, 'H': 1, 'R': 1, 'S': 1, 'T': 1, 'D': 1})
        self.assertEqual(Ficha4.contar_bases('8172bfw9df8qwolwof997'), 'Introduza uma sequência de DNA, RNA ou de aminoácidos válida')
        
    def reading_frames(self):
        self.assertEqual(Ficha4.reading_frames('ATGCATAGGTCGACCGATTGA'), ['MHRSTD_', 'CIGRPI', 'A_VDRL', 'SIGRPMH', 'QSVDLC', 'NRSTYA']
        self.assertEqual(Ficha4.reading_frames('atgcataggtcgaccgattga'), ['MHRSTD_', 'CIGRPI', 'A_VDRL', 'SIGRPMH', 'QSVDLC', 'NRSTYA']
        self.assertEqual(Ficha4.reading_frames('snd1268371iehefowdq91'), 'Introduza uma sequência de DNA válida')

    def get_proteins(self):
        self.assertEqual(Ficha4.get_proteins('ATGGCAGATTCCAACGGTACTATTACCGTTGAAGAGCTTAAAAAGCTCCTTGAACAATGGAACCTAGTAATAGGTTTCCTATTCCTTACATGGATTTGTCTTCTACAATTTGCCTATGCCAACAGGAATAGGTTTTTGTATATAATTAAGTTAATTTTCCTCTGGCTGTTATGGCCAGTAACTTTAGCTTGTTTTGTGCTTGCTGCTGTTTACAGAATAAATTGGATCACCGGTGGAATTGCTATCGCAATGGCTTGTCTTGTAGGCTTGATGTGGCTCAGCTACTTCATTGCTTCTTTCAGACTGTTTGCGCGTACGCGTTCCATGTGGTCATTCAATCCAGAAACTAACATTCTTCTCAACGTGCCACTCCATGGCACTATTCTGACCAGACCGCTTCTAGAAAGTGAACTCGTAATCGGAGCTGTGATCCTTCGTGGACATCTTCGTATTGCTGGACACCATCTAGGACGCTGTGACATCAAGGACCTGCCTAAAGAAATCACTGTTGCTACATCACGAACGCTTTCTTATTACAAATTGGGAGCTTCGCAGCGTGTAGCAGGTGACTCAGGTTTTGCTGCATACAGTCGCTACAGGATTGGCAACTATAAATTAAACACAGACCATTCCAGTAGCAGTGACAATATTGCTTTGCTTGTACAGTAA'), ['MADSNGTITVEELKKLLEQWNLVIGFLFLTWICLLQFAYANRNRFLYIIKLIFLWLLWPVTLACFVLAAVYRINWITGGIAIAMACLVGLMWLSYFIASFRLFARTRSMWSFNPETNILLNVPLHGTILTRPLLESELVIGAVILRGHLRIAGHHLGRCDIKDLPKEITVATSRTLSYYKLGASQRVAGDSGFAAYSRYRIGNYKLNTDHSSSSDNIALLVQ_', 'MSQRPRWCPAIRRCPRRITAPITSSLSRSGLVRIVPWSGTLRRMLVSGLNDHMERVRANSLKEAMK_', 'MASNFSLFCACCCLQNKLDHRWNCYRNGLSCRLDVAQLLHCFFQTVCAYAFHVVIQSRN_', 'MEPSNRFPIPYMDLSSTICLCQQE_', 'MQQNLSHLLHAAKLPICNKKAFVM_', 'MVSSNTKMSTKDHSSDYEFTF_', 'MEWHVEKNVSFWIE_', 'MTTWNAYAQTV_', 'MPTGIGFCI_', 'MALF_', 'MVCV_', 'M_'])
        self.assertEqual(Ficha4.get_proteins('sljfbwiyf3746192yiofncq2'), 'Introuza uma sequência de DNA válida')
                         
if __name__ == '__main__':
    unittest.main()
