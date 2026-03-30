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
int dataset_size = 0;

void show_menu() {
    printf("\n==== DNA Matching System ====\n");
    printf("1. Load Sample Data\n");
    printf("2. Enter Query DNA Sequence\n");
    printf("3. Validate Query\n");
    printf("4. Exit\n");
    printf("5. Search Exact Match (Trie)\n");
    printf("6. Search Pattern (Suffix Tree)\n");
    printf("7. Rank Matches (Skip List)\n");
    printf("Choose an option: ");
}

// Load file and store sequences
void load_and_store(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Error opening %s\n", filename);
        return;
    }

    char line[MAX_SEQ];

    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = 0;

        // Skip labels
        if (line[0] == '>') continue;

        normalize_sequence(line);

        if (!validate_sequence(line)) continue;

        // Store in dataset
        strcpy(dataset[dataset_size], line);
        dataset_size++;

        // Insert into Trie
        insert_sequence(trie_root, line);
    }

    fclose(file);
}

int calculate_similarity(const char* seq, const char* query) {
    
    int match = 0;
    int len1 = strlen(seq);
    int len2 = strlen(query);
    int min = len1 < len2 ? len1 : len2;
    if (min == 0) return 0;

    for (int i = 0; i < min; i++) {
        if (seq[i] == query[i]) {
            match++;
        }
    }
    return (match * 100) / min; // percentage
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

                load_and_store("data/human.txt");
                load_and_store("data/chimpanzee.txt");
                load_and_store("data/mouse.txt");
                load_and_store("data/virus_strain_a.txt");
                load_and_store("data/virus_strain_b.txt");

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
                        printf("Pattern found in sequence %d\n", i + 1);
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
                    printf("Load dataset first.\n");
                    break;
                }

                if (!validate_sequence(query)) {
                    printf("Enter a valid query first.\n");
                    break;
                }
                
                // Clear previous skip list
                free_skiplist();
                init_skiplist();   // 🔥 MUST
                printf("Calculating similarity scores...\n");

                for (int i = 0; i < dataset_size; i++) {
                    if (strlen(dataset[i]) == 0) continue;
                    int score = calculate_similarity(dataset[i], query);
                    printf("Seq %d Score: %d\n", i, score);  // debug
                    if (score > 0) {
                        insert_skiplist(dataset[i], score);
                    }
                }
                display_top_matches(5); // Top 5 matches
                break;
            }
            default:
                printf("Invalid choice.\n");
        }
    }
    return 0;
}
