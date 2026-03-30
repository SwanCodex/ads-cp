#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "trie.h"
#include "suffix_tree.h"
#include "skiplist.h"

#define MAX_SEQ 1000
#define MAX_DATASET 100

TrieNode* trie_root = NULL;

// Store dataset in memory
char dataset[MAX_DATASET][MAX_SEQ];
char species[MAX_DATASET][50];
char labels[MAX_DATASET][50];
int dataset_size = 0;

void show_menu() {
    printf("\n========================================\n");
    printf("        🧬 DNA MATCHING SYSTEM\n");
    printf("========================================\n");
    printf("1. Load DNA Dataset\n");
    printf("2. Enter DNA Query\n");
    printf("3. Validate Query\n");
    printf("4. Exit\n");
    printf("5. Exact Match (Trie)\n");
    printf("6. Pattern Search (Suffix Tree)\n");
    printf("7. Rank & Analyze Matches\n");
    printf("========================================\n");
    printf("Enter choice: ");
}

// Load file and store sequences
void load_and_store(const char* filename, const char* species_name) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening %s\n", filename);
        return;
    }

    char line[MAX_SEQ];
    char current_label[50];

    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;

        if (line[0] == '>') {
            strcpy(current_label, line + 1); // remove '>'
            continue;
        }

        normalize_sequence(line);
        if (!validate_sequence(line)) continue;

        strcpy(dataset[dataset_size], line);
        strcpy(species[dataset_size], species_name);
        strcpy(labels[dataset_size], current_label);

        insert_sequence(trie_root, line);

        dataset_size++;
    }

    fclose(file);
}

int calculate_similarity(const char* seq, const char* query) {
    int len1 = strlen(seq);
    int len2 = strlen(query);

    int max_match = 0;

    // Slide query across sequence
    for (int i = 0; i <= len1 - len2; i++) {
        int match = 0;

        for (int j = 0; j < len2; j++) {
            if (seq[i + j] == query[j]) {
                match++;
            }
        }

        if (match > max_match) {
            max_match = match;
        }
    }

    return (max_match * 100) / len2;
}

int main() {
    int choice;
    char query[MAX_SEQ];

    while (1) {
        show_menu();
        scanf("%d", &choice);
        getchar();

        switch (choice) {

            case 1:
                printf("\nLoading dataset...\n");

                // Reset
                dataset_size = 0;

                if (trie_root) free_trie(trie_root);
                trie_root = create_trie();

                load_and_store("data/human.txt", "Human");
                load_and_store("data/chimpanzee.txt", "Chimpanzee");
                load_and_store("data/mouse.txt", "Mouse");
                load_and_store("data/virus_strain_a.txt", "Virus A");
                load_and_store("data/virus_strain_b.txt", "Virus B");

                printf("Dataset loaded: %d sequences\n", dataset_size);
                break;

            case 2:
                printf("\nEnter DNA query: ");
                fgets(query, MAX_SEQ, stdin);
                query[strcspn(query, "\n")] = 0;

                normalize_sequence(query);
                printf("Normalized Query: %s\n", query);
                break;

            case 3:
                if (validate_sequence(query)) {
                    printf("Valid DNA sequence.\n");
                } else {
                    printf("Invalid DNA sequence.\n");
                }
                break;

            case 4:
                printf("Exiting...\n");
                return 0;

            case 5:
                if (!trie_root) {
                    printf("Load dataset first.\n");
                    break;
                }

                if (search_sequence(trie_root, query)) {
                    printf("Exact sequence FOUND in dataset.\n");
                } else {
                    printf("Exact sequence NOT found.\n");
                }
                break;

            case 6: {
                if (dataset_size == 0) {
                    printf("Load dataset first.\n");
                    break;
                }

                int found = 0;

                printf("Searching pattern across dataset...\n");

                for (int i = 0; i < dataset_size; i++) {
                    build_suffix_tree(dataset[i]);

                    if (search_pattern(query)) {
                        printf("✔ Pattern found in: %s (%s)\n",
                            species[i],
                            labels[i]);
                        found = 1;
                    }

                    free_suffix_tree();
                }

                if (!found) {
                    printf("Pattern NOT found in dataset.\n");
                }

                break;
            }
            case 7: {
                if (dataset_size == 0) {
                    printf("⚠️ Load dataset first.\n");
                    break;
                }

                if (!validate_sequence(query)) {
                    printf("⚠️ Enter a valid query first.\n");
                    break;
                }

                free_skiplist();
                init_skiplist();

                printf("\n🔬 Analyzing DNA Sequence...\n");
                printf("----------------------------------------\n");

                int best_score = 0;
                char best_species[50];
                strcpy(best_species, "None");

                for (int i = 0; i < dataset_size; i++) {

                    int score = calculate_similarity(dataset[i], query);

                    if (score > 0) {
                        insert_skiplist(dataset[i], species[i], score);                    
                    }

                    if (score > best_score) {
                        best_score = score;
                        strcpy(best_species, species[i]);
                    }
                }

                printf("\n🧬 BEST MATCH\n");
                printf("----------------------------------------\n");
                printf("Species: %s\n", best_species);
                printf("Similarity: %d%%\n", best_score);

                printf("\n🏆 TOP MATCHES\n");
                printf("----------------------------------------\n");
                display_top_matches(5);

                break;
            }
            default:
                printf("Invalid choice.\n");
        }
    }

    return 0;
}