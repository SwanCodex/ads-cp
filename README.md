# 🧬 DNA Analysis System (C Project)

## 🚀 How to Run

```bash
gcc *.c -o dna_app
./dna_app
```

## 📌 Overview

This project is a **DNA analysis system** focused on real-world biological insights rather than simple sequence matching.

It answers key questions such as:
- Which species does this DNA belong to?
- What mutations are present?
- Is the sequence viral or normal?
- What is the closest known gene?

---

## 🧠 Core Features

### 🔍 1. Species Identification
- Determines most likely organism (Human / Virus / Mouse / Chimpanzee)
- Based on global alignment scoring

---

### 🧬 2. Mutation Detection
- Detects base-level mutations
- Reports:
  - Position
  - Original base → mutated base
  - Total mutation count

---

### ⚠️ 3. Disease / Viral Detection
- Identifies whether sequence resembles viral DNA
- Flags potential abnormal sequences

---

### 📊 4. Top Match Ranking
- Uses Skip List to rank best matches
- Displays top 3 closest sequences

---

### 🧪 5. Alignment Visualization
- Uses Needleman-Wunsch algorithm
- Shows:
  - Matches (|)
  - Mismatches (*)
  - Gaps

---

### ⚡ 6. Fast Filtering (K-mer Hashing)
- Eliminates unrelated sequences early
- Improves performance significantly

---

## 🧱 Data Structures Used

| Structure        | Role |
|----------------|------|
| Trie           | Stores DNA sequences |
| Hash Table     | K-mer filtering |
| Suffix Tree    | Pattern detection |
| Skip List      | Ranking results |
| DP Matrix      | Alignment (Needleman-Wunsch) |

---

## 🔬 Workflow

1. Load dataset (multi-species DNA)
2. Input query DNA
3. Filter using k-mers
4. Perform alignment
5. Identify species
6. Detect mutations
7. Rank top matches
8. Display alignment

---