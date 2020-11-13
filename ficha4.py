def ler_seq(FileHandle):
    """

Devolve uma sequência de nucleótidos a partir de um ficheiro aberto
    
PARAMETERS
----------
FiheHandle: 
return: COnteúdo do ficheiro
    """
    ficheiro = open('sequence.txt', 'r')
    conteudo = ficheiro.read()
    return conteudo

def ler_FASTA_seq(FileHandle):
    """

Devolve uma sequência a partir de um ficheiro aberto

PARAMETERS
----------
FileHandle: 
return: Sequência presente no ficheiro sob formato FASTA
    """
    with open('sequence.fasta') as FH:
        while (s := ler_FASTA_seq(FH)):
            print(s)

def complemento_inverso(seq):
    """

Devolve a sequência complementar e inversa de uma sequência de DNA

PARAMETERS
----------
seq: str
    Uma sequência de nucleótidos constituida apenas por A, T, G e C,
    tanto em caracteres maiúsculos como minúsculos
return: Uma sequência complementar à primeira (A->T, T->A, G->C e C->G),
invertendo depois essa que seria 3'->5' para 5'->3', em caracteteres maiúsculos.
    """
    seq = seq.lower()
    seq1 = seq[::-1]
    seq2 = seq1.replace('a','T').replace('t','A').replace('g','C').replace('c','G')
    return seq2.upper()

def transcricao(seq):
    """

Efetua a transcrição de uma sequência de DNA

PARAMETERS
----------
seq: str
    Uma sequência de nucleótidos constituida apenas por A, T, G e C,
    tanto em caracteres maiúsculos como minúsculos
return: Uma string com a sequência transcrita em maiúsculas
    """
    seq = seq.upper()
    seq1= seq.replace('T','U')
    return seq1

def traducao(seq):
    """

Efetua a tradução de uma sequência de DNA

PARAMETERS
----------
seq: str
    Uma sequência de DNA
return: Uma string com a sequência de aminoácidos correspondentes à sequência introduzida

    """
    gencode = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    import re
    seq = seq.upper()
    codao = re.findall('...',seq)
    amino = []
    for x in codao:
        amino.append(gencode[x])
    resultado = ''.join(amino)
    return resultado

def valida(seq):
    """

Diz se o input é considerado uma sequência de DNA ou não

PARAMETERS
----------
seq: str
    Uma sequência de caracteres
return: TRUE se a sequência providenciada corresponder a uma sequência de DNA
com A, T, C e G; FALSE caso contrário
    """
    seq = seq.upper()
    if len(set(seq) - {'A','T','C','G'}) == 0:
        return "TRUE"
    else:
        return "FALSE"
    
def contar_bases(seq):
    """

Devolve o número de vezes que cada nucleótido aparece na sequência de DNA/RNA ou aminoácido no caso se ser uma sequência proteica
    
PARAMETERS
----------

seq: str
   Uma sequência de DNA, RNA ou aminiácidos     
return: Dicionário com a frequência de cada nucleótido/aminoácido presente na sequência
    """
    seq = seq.upper()
    letras={}
    for letra in sorted(seq):
        if letra not in letras:
            letras[letra] = 0
        letras[letra] += 1
    return letras

def reading_frames(seq):
    """

Devolve as reading frames dentro da sequência de DNA

PARAMETERS
----------

seq: str
    Uma sequência de DNA
return: Uma lista com 3 reading frames de 5'->3' e mais 3 com o complemento inverso 3'->5' da cadeia de DNA
    """
    rf = []
    for y in range(3):
        rf.append(''.join(traducao(seq[y:])))
    for y in range(3):
        rf.append(''.join(traducao(complemento_inverso(seq)[y:])))
    return rf
    
def get_proteins(seq):
    """

Devolve uma lista com todas as proteínas possiveis dentro da sequência de DNA
    
PARAMETERS
----------

seq: str
    Uma sequência de DNA
return: Uma lista com todas as proteinas encontradas nas reading frames da sequência de DNA,
ordenadas por ordem decrescente de tamanho e dentro do mesmo tamanho, por ordem alfabética
    """
    import re
    proteinas = [re.findall("M.*?_", orfs) for orfs in reading_frames(seq)]
    proteinassorted = sorted({p for lp in proteinas for p in lp}, key = lambda x: (-len(x), x))
    return proteinassorted
