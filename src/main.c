#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "utils.h"
#include "trie.h"
#include "suffix_tree.h"
#include "skiplist.h"
#include "alignment.h"
#include "hashtable.h"

#define MAX_SEQ 2000
#define MAX_DATASET 200
#define MAX_DISEASE 100
#define KMER 4

TrieNode* trie_root = NULL;

/* ================= DATA ================= */
char dataset[MAX_DATASET][MAX_SEQ];
char species[MAX_DATASET][50];
char labels[MAX_DATASET][50];
int dataset_size = 0;

/* ================= DISEASE DB ================= */
char disease_seq[MAX_DISEASE][MAX_SEQ];
char disease_name[MAX_DISEASE][50];
int disease_count = 0;

/* ================= LOAD DATASET ================= */
void load_and_store(const char* filename, const char* species_name) {
    FILE* file = fopen(filename, "r");
    if (!file) return;

    char line[MAX_SEQ];
    char current_label[50] = "";

    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;

        if (strlen(line) == 0) continue;

        if (line[0] == '>') {
            strcpy(current_label, line + 1);
            continue;
        }

        normalize_sequence(line);
        if (!validate_sequence(line)) continue;

        strcpy(dataset[dataset_size], line);
        strcpy(species[dataset_size], species_name);
        strcpy(labels[dataset_size], current_label);

        insert_sequence(trie_root, line);

        int len = strlen(line);
        for (int i = 0; i <= len - KMER; i++) {
            char kmer[KMER + 1];
            strncpy(kmer, &line[i], KMER);
            kmer[KMER] = '\0';
            insert_kmer(kmer);
        }

        dataset_size++;
    }

    fclose(file);
}

/* ================= LOAD DISEASE DB ================= */
void load_disease_db(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) return;

    char line[MAX_SEQ];
    char current_label[50] = "";

    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;

        if (strlen(line) == 0) continue;

        if (line[0] == '>') {
            strcpy(current_label, line + 1);
            continue;
        }

        normalize_sequence(line);
        if (!validate_sequence(line)) continue;

        strcpy(disease_seq[disease_count], line);
        strcpy(disease_name[disease_count], current_label);
        disease_count++;
    }

    fclose(file);
}

/* ================= SUFFIX TREE HELPER ================= */
int suffix_tree_match(const char* sequence, const char* query) {
    if (strlen(sequence) == 0 || strlen(query) == 0)
        return 0;

    build_suffix_tree(sequence);
    int found = search_pattern(query);
    free_suffix_tree();

    return found;
}

/* ================= MUTATION DETECTION ================= */
void detect_mutations(const char* ref, const char* query) {
    if (strlen(ref) == 0 || strlen(query) == 0) return;

    int len = strlen(ref) < strlen(query) ? strlen(ref) : strlen(query);
    int mutations = 0;

    printf("\nMutation Analysis:\n");

    for (int i = 0; i < len; i++) {
        if (ref[i] != query[i]) {
            printf("Position %d: %c -> %c\n", i + 1, ref[i], query[i]);
            mutations++;
        }
    }

    printf("Total Mutations: %d\n", mutations);
}

/* ================= SPECIES IDENTIFICATION ================= */
int find_best_species(const char* query) {
    int best_score = -999999;
    int best_index = -1;

    for (int i = 0; i < dataset_size; i++) {
        if (strlen(dataset[i]) == 0) continue;

        int score = needleman_wunsch(dataset[i], query);

        if (suffix_tree_match(dataset[i], query)) {
            score += 20;
        }

        if (score > best_score) {
            best_score = score;
            best_index = i;
        }
    }

    if (best_index != -1) {
        printf("\nMost Likely Origin:\n");
        printf("Species: %s\n", species[best_index]);
        printf("Gene   : %s\n", labels[best_index]);
        printf("Score  : %d\n", best_score);

        if (strstr(species[best_index], "Virus"))
            printf("Possible viral sequence detected.\n");
        else
            printf("Likely normal organism DNA.\n");
    }

    return best_index;
}

/* ================= DISEASE ANALYSIS ================= */
void analyze_disease_risk(const char* query) {
    double best_score = -9999;
    int best_index = -1;

    for (int i = 0; i < disease_count; i++) {
        if (strlen(disease_seq[i]) == 0) continue;

        int raw_score = needleman_wunsch(disease_seq[i], query);

        int len1 = strlen(disease_seq[i]);
        int len2 = strlen(query);

        double norm_score = (double)raw_score / ((len1 + len2) / 2.0);

        if (norm_score > best_score) {
            best_score = norm_score;
            best_index = i;
        }
    }

    if (best_index == -1) return;

    printf("\n=== Genetic Risk Analysis ===\n");
    printf("Closest Marker: %s\n", disease_name[best_index]);

    int risk = (int)((best_score + 2) / 3 * 100);

    if (risk < 0) risk = 0;
    if (risk > 100) risk = 100;

    printf("Risk Score: %d%%\n", risk);

    if (risk > 85)
        printf("HIGH RISK\n");
    else if (risk > 60)
        printf("MODERATE RISK\n");
    else
        printf("LOW RISK\n");
}

/* ================= MAIN ================= */
int main() {
    char query[MAX_SEQ];

    printf("\n=== DNA ANALYSIS SYSTEM ===\n");

    trie_root = create_trie();
    free_table();

    load_and_store("data/human.txt", "Human");
    load_and_store("data/chimpanzee.txt", "Chimpanzee");
    load_and_store("data/mouse.txt", "Mouse");
    load_and_store("data/virus_strain_a.txt", "Virus A");
    load_and_store("data/virus_strain_b.txt", "Virus B");

    load_disease_db("data/disease_db.txt");

    printf("Dataset loaded: %d sequences\n", dataset_size);

    /* INPUT */
    while (1) {
        printf("\nEnter DNA sequence (-1 to exit): ");
        fgets(query, MAX_SEQ, stdin);
        if(query[0] == '-'){
            break;
        }
        query[strcspn(query, "\n")] = 0;

        normalize_sequence(query);

        if (strlen(query) == 0) {
            printf("Empty input not allowed.\n");
            continue;
        }

        if (!validate_sequence(query)) {
            printf("Invalid DNA sequence.\n");
            continue;
        }

    /* HASH FILTER */
    if (strlen(query) >= KMER) {
        char kmer[KMER + 1];
        strncpy(kmer, query, KMER);
        kmer[KMER] = '\0';

        if (!search_kmer(kmer)) {
            printf("\nNo biological similarity found.\n");
            return 0;
        }
    }

    printf("\n=== ANALYSIS RESULTS ===\n");

    if (search_sequence(trie_root,query)){
        printf("exact sequence found!!");
    }
    else{
    /* SUFFIX TREE PATTERN DETECTION */
    printf("\nPattern Detection (Suffix Tree):\n");
    int any_pattern = 0;

    for (int i = 0; i < dataset_size; i++) {
        if (suffix_tree_match(dataset[i], query)) {
            printf("Pattern found in: %s (%s)\n", species[i], labels[i]);
            any_pattern = 1;
        }
    }

    if (!any_pattern) {
        printf("No pattern found.\n");
    }
}

    int best_index = find_best_species(query);

    /* RANKING */
    init_skiplist();
    for (int i = 0; i < dataset_size; i++) {
        if (strlen(dataset[i]) == 0) continue;

        int score = needleman_wunsch(dataset[i], query);

        if (suffix_tree_match(dataset[i], query)) {
            score += 20;
        }

        insert_skiplist(dataset[i], species[i], score);
    }

    printf("\nTop Matches:\n");
    display_top_matches(3);

    /* MUTATION */
    if (best_index != -1)
        detect_mutations(dataset[best_index], query);

    /* DISEASE */
    analyze_disease_risk(query);

    /* ALIGNMENT */
    if (best_index != -1)
        print_alignment(dataset[best_index], query);
}


    printf("\n=== END ===\n");


    return 0;
}