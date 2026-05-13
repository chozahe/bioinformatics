#!/usr/bin/env python3
"""
Задание №1: Поиск длиннейшей ORF (Open Reading Frame)
Генерация случайной ДНК с заданным GC-составом и поиск самой длинной ORF
во всех 6 рамках трансляции.
"""

import random
import sys

# ─── Генетический код (однобуквенные обозначения) ───
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# Трёхбуквенные обозначения аминокислот
THREE_LETTER = {
    'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu',
    'F': 'Phe', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn',
    'P': 'Pro', 'Q': 'Gln', 'R': 'Arg', 'S': 'Ser',
    'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr',
}

START_CODON = 'ATG'
STOP_CODONS = {'TAA', 'TAG', 'TGA'}


# ─── Вспомогательные функции ───

def input_int(prompt, min_val, max_val):
    """Безопасный ввод целого числа с проверкой диапазона."""
    while True:
        try:
            val = int(input(prompt))
            if min_val <= val <= max_val:
                return val
            print(f"Ошибка: введите число от {min_val} до {max_val}.")
        except ValueError:
            print("Ошибка: введите целое число.")


def generate_dna(length, gc_content):
    """
    Генерирует случайную строку ДНК заданной длины и GC-состава.
    GC-состав задаётся в процентах (20-80). 
    G% = GC%/2, C% = GC%/2, A% = T% = (100-GC%)/2.
    """
    g_percent = gc_content / 2.0
    c_percent = gc_content / 2.0
    a_percent = (100 - gc_content) / 2.0
    t_percent = (100 - gc_content) / 2.0

    nucleotides = ['A', 'T', 'G', 'C']
    weights = [a_percent, t_percent, g_percent, c_percent]

    return ''.join(random.choices(nucleotides, weights=weights, k=length))


def reverse_complement(dna):
    """Возвращает обратно-комплементарную цепь ДНК."""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(dna))


def translate_codon(codon):
    """Переводит кодон в аминокислоту (однобуквенный код)."""
    return GENETIC_CODE.get(codon.upper(), 'X')


def translate_sequence(seq, one_letter=True):
    """Переводит последовательность нуклеотидов в белок."""
    proteins = []
    for i in range(0, len(seq), 3):
        if i + 3 <= len(seq):
            aa = translate_codon(seq[i:i+3])
            if one_letter:
                proteins.append(aa)
            else:
                proteins.append(THREE_LETTER.get(aa, aa))
    return ' '.join(proteins) if one_letter else ' '.join(proteins)


def format_triplets(seq):
    """Разделяет последовательность на триплеты через пробел."""
    return ' '.join(seq[i:i+3] for i in range(0, len(seq), 3))


# ─── Поиск ORF ───

def find_longest_orf_in_sequence(sequence, frame_offset):
    """
    Ищет самую длинную ORF в заданной последовательности 
    в указанной рамке считывания (0, 1 или 2).
    
    Возвращает (start_pos, end_pos, orf_sequence) или None.
    Позиции — 1-индексированные относительно переданной последовательности.
    """
    n = len(sequence)
    longest_orf = None  # (start, end, seq)

    # Сканируем последовательность, начиная со смещения frame_offset
    i = frame_offset
    while i + 3 <= n:
        codon = sequence[i:i+3]

        # Ищем старт-кодон
        if codon == START_CODON:
            orf_start = i
            # Идём от старт-кодона до стоп-кодона шагами по 3
            j = i
            while j + 3 <= n:
                c = sequence[j:j+3]
                if c in STOP_CODONS:
                    orf_end = j + 3  # включительно
                    orf_seq = sequence[orf_start:orf_end]
                    # Проверяем минимальную длину белка (>= 10 аминокислот)
                    protein_length = len(orf_seq) // 3 - 1  # -1 за стоп-кодон
                    if protein_length >= 10:
                        if longest_orf is None or len(orf_seq) > len(longest_orf[2]):
                            longest_orf = (orf_start, orf_end, orf_seq)
                    break  # выходим, даже если условие не выполнено
                j += 3
            # Перемещаем i на позицию после старт-кодона,
            # чтобы найти другие ORF
            i = orf_start + 3
        else:
            i += 3

    return longest_orf


def find_longest_orf(sequence, is_reverse=False):
    """
    Ищет самую длинную ORF во всех 3 рамках считывания последовательности.
    
    Возвращает (orf_seq, start_pos, end_pos, frame, chain) или None.
    start_pos, end_pos — позиции относительно всей исходной ДНК (1-индексированные).
    """
    chain_name = "обратная (обратно-комплементарная)" if is_reverse else "прямая"
    best_orf = None

    for frame in range(3):
        result = find_longest_orf_in_sequence(sequence, frame)
        if result is not None:
            start, end, orf_seq = result
            # Позиции: 1-индексированные в рамках этой последовательности
            # start - это индекс в sequence (0-based), frame - смещение
            # Реальная позиция в sequence: start+1 (1-индексированная)
            abs_start = start + 1
            abs_end = end  # end уже exclusive-индекс, но мы хотим включительно
            if best_orf is None or len(orf_seq) > len(best_orf[0]):
                best_orf = (orf_seq, abs_start, abs_end, frame + 1, chain_name)

    return best_orf


def find_longest_orf_all_frames(dna_sequence):
    """
    Ищет самую длинную ORF во всех 6 рамках трансляции.
    
    Возвращает результат в виде словаря или None.
    """
    # Прямая цепь — 3 рамки
    fwd_result = find_longest_orf(dna_sequence, is_reverse=False)
    
    # Обратно-комплементарная цепь — 3 рамки
    rev_seq = reverse_complement(dna_sequence)
    rev_result = find_longest_orf(rev_seq, is_reverse=True)

    # Сравниваем результаты
    best = None
    
    if fwd_result is not None:
        best = {
            'orf_seq': fwd_result[0],
            'start': fwd_result[1],
            'end': fwd_result[2],
            'frame': fwd_result[3],
            'chain': fwd_result[4],
            'is_reverse': False,
        }
    
    if rev_result is not None:
        rev_best = {
            'orf_seq': rev_result[0],
            'start': rev_result[1],
            'end': rev_result[2],
            'frame': rev_result[3],
            'chain': rev_result[4],
            'is_reverse': True,
        }
        if best is None or len(rev_best['orf_seq']) > len(best['orf_seq']):
            best = rev_best

    return best


def get_orf_protein_info(orf_seq):
    """Извлекает информацию о белке из ORF."""
    # Отрезаем стоп-кодон для длины белка
    coding_seq = orf_seq[:-3]  # без стоп-кодона
    orf_triplets = format_triplets(orf_seq)
    
    # Однобуквенные обозначения
    protein_one = []
    for i in range(0, len(coding_seq), 3):
        codon = coding_seq[i:i+3]
        protein_one.append(translate_codon(codon))
    
    # Трёхбуквенные обозначения
    protein_three = []
    for aa in protein_one:
        protein_three.append(THREE_LETTER.get(aa, aa))
    
    return {
        'orf_triplets': orf_triplets,
        'protein_one_letter': ' '.join(protein_one),
        'protein_three_letter': ' '.join(protein_three),
        'protein_length': len(protein_one),
    }


# ─── Основной блок ───

def main():
    print("=" * 60)
    print("   Задание №1: Поиск длиннейшей ORF")
    print("=" * 60)
    print()

    # Ввод параметров
    print("Введите параметры ДНК:")
    length = input_int("  Длина (100-1000): ", 100, 1000)
    gc = input_int("  GC-состав в % (20-80): ", 20, 80)

    print()
    print(f"Генерирую ДНК: длина={length}, GC={gc}%...")
    
    # Генерация ДНК
    random.seed()  # инициализация генератора
    dna = generate_dna(length, gc)

    print()
    print("─" * 60)
    print("Сгенерированная ДНК (прямая цепь):")
    print("─" * 60)
    print(dna)
    
    # GC-состав сгенерированной строки
    actual_gc = (dna.count('G') + dna.count('C')) / len(dna) * 100
    print(f"\nФактический GC-состав: {actual_gc:.1f}%")
    print(f"Содержание: A={dna.count('A')} ({dna.count('A')/len(dna)*100:.1f}%), "
          f"T={dna.count('T')} ({dna.count('T')/len(dna)*100:.1f}%), "
          f"G={dna.count('G')} ({dna.count('G')/len(dna)*100:.1f}%), "
          f"C={dna.count('C')} ({dna.count('C')/len(dna)*100:.1f}%)")

    # Поиск ORF
    print()
    print("─" * 60)
    print("Поиск длиннейшей ORF во всех 6 рамках трансляции...")
    print("─" * 60)
    
    result = find_longest_orf_all_frames(dna)

    if result is None:
        print()
        print("❌ ORF не найдена: ни одна ORF с длиной белка ≥ 10 а.о. не обнаружена.")
        return

    orf_seq = result['orf_seq']
    protein_info = get_orf_protein_info(orf_seq)

    print()
    print("✅ ДЛИННЕЙШАЯ ORF НАЙДЕНА!")
    print()

    # Если ORF на обратной цепи, пересчитываем координаты
    if result['is_reverse']:
        # result['start'] и result['end'] — позиции на обратной цепи
        # Пересчитываем на прямую цепь: 
        # position_on_forward = length - position_on_reverse + 1
        rev_start_on_fwd = length - result['end'] + 1
        rev_end_on_fwd = length - result['start'] + 1
        
        print(f"  Цепь: {result['chain']}")
        print(f"  Рамка считывания: {result['frame']}")
        print(f"  Позиции на обратной цепи: {result['start']}–{result['end']}")
        print(f"  Позиции на прямой цепи (reverse complement): {rev_start_on_fwd}–{rev_end_on_fwd}")
    else:
        print(f"  Цепь: {result['chain']}")
        print(f"  Рамка считывания: {result['frame']}")
        print(f"  Позиции в цепи: {result['start']}–{result['end']}")
    
    print(f"  Длина ORF: {len(orf_seq)} нуклеотидов")
    print(f"  Длина белка: {protein_info['protein_length']} аминокислот")
    
    print()
    print("─" * 60)
    print("ORF с разделением на триплеты:")
    print("─" * 60)
    print(protein_info['orf_triplets'])
    
    print()
    print("─" * 60)
    print("Транслированный белок (однобуквенные обозначения):")
    print("─" * 60)
    print(protein_info['protein_one_letter'].upper())
    
    print()
    print("─" * 60)
    print("Транслированный белок (трёхбуквенные обозначения):")
    print("─" * 60)
    print(protein_info['protein_three_letter'])
    print()

    # Вывод всей ДНК с подсвеченной ORF (для проверки)
    print("─" * 60)
    print("ДНК с найденной ORF (в квадратных скобках):")
    print("─" * 60)
    
    if result['is_reverse']:
        rev_comp = reverse_complement(dna)
        start = result['start'] - 1
        end = result['end']
        annotated = rev_comp[:start] + "[" + rev_comp[start:end] + "]" + rev_comp[end:]
        print(annotated)
    else:
        start = result['start'] - 1
        end = result['end']
        annotated = dna[:start] + "[" + dna[start:end] + "]" + dna[end:]
        print(annotated)

    print()
    print("=" * 60)
    print("   Готово! Рекомендуется проверить результат")
    print("   в ORFfinder (https://www.ncbi.nlm.nih.gov/orffinder/)")
    print("=" * 60)


if __name__ == '__main__':
    main()
