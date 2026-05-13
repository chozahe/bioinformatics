#!/usr/bin/env python3
"""
Задание №2: Сбор статистик FASTQ-файла
Принимает FASTQ-файл (Phred33), выводит статистики.
"""

import sys
import gzip
import os


def open_fastq(filename):
    """Открывает FASTQ-файл (обычный или gzip)."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')


def count_lines(filename):
    """Быстрый подсчёт строк (для больших файлов)."""
    with open_fastq(filename) as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


def parse_fastq(filename):
    """
    Генератор, читает FASTQ и выдаёт (header, sequence, quality) для каждого рида.
    """
    with open_fastq(filename) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            sep = f.readline().strip()
            qual = f.readline().strip()
            
            if not header.startswith('@'):
                # На случай если файл кривой — пытаемся восстановиться
                continue
            
            yield header, seq, qual


def compute_stats(filename):
    """
    Собирает статистики FASTQ-файла.
    """
    read_count = 0
    total_bases = 0
    total_quality_bases = 0  # bases with Q >= 30
    
    for header, seq, qual in parse_fastq(filename):
        read_count += 1
        seq_len = len(seq)
        total_bases += seq_len
        
        # Phred33: Q = ord(symbol) - 33
        # Q >= 30 → ord(symbol) >= 63 (т.е. символ ASCII >= '?')
        for ch in qual:
            if ord(ch) >= 63:
                total_quality_bases += 1
    
    return read_count, total_bases, total_quality_bases


def main():
    if len(sys.argv) < 2:
        print("Использование: python3 fastq_stats.py <filename.fastq[.gz]>")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    if not os.path.exists(filename):
        print(f"Ошибка: файл '{filename}' не найден.")
        sys.exit(1)
    
    print(f"Обрабатываю: {filename}")
    print()
    
    read_count, total_bases, total_quality_bases = compute_stats(filename)
    
    avg_length = total_bases / read_count if read_count > 0 else 0
    quality_percent = (total_quality_bases / total_bases * 100) if total_bases > 0 else 0
    
    print("═" * 40)
    print("     СТАТИСТИКИ FASTQ-ФАЙЛА")
    print("═" * 40)
    print(f"  1. Количество ридов:          {read_count}")
    print(f"  2. Общее количество букв:     {total_bases}")
    print(f"  3. Средняя длина рида:        {avg_length:.1f}")
    print(f"  4. Процент Q ≥ 30:            {quality_percent:.2f}%")
    print("═" * 40)


if __name__ == '__main__':
    main()
