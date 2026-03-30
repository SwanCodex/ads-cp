#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hashtable.h"

// Simple hash function (djb2)
unsigned long hash(const char* str) {
    unsigned long hash = 5381;
    int c;
    while ((c = *str++)) {
        hash = ((hash << 5) + hash) + c; // hash * 33 + c
    }
    return hash % TABLE_SIZE;
}

static HashNode* table[TABLE_SIZE] = {NULL};

void insert_kmer(const char* kmer) {
    unsigned long index = hash(kmer);
    
    // Check if duplicate (optional but helpful for k-mers)
    HashNode* current = table[index];
    while (current) {
        if (strcmp(current->kmer, kmer) == 0) return;
        current = current->next;
    }

    // Insert at head of list (chaining)
    HashNode* newNode = (HashNode*)malloc(sizeof(HashNode));
    newNode->kmer = strdup(kmer);
    newNode->next = table[index];
    table[index] = newNode;
}

int search_kmer(const char* kmer) {
    unsigned long index = hash(kmer);
    HashNode* current = table[index];
    while (current) {
        if (strcmp(current->kmer, kmer) == 0) return 1;
        current = current->next;
    }
    return 0;
}

void free_table() {
    for (int i = 0; i < TABLE_SIZE; i++) {
        HashNode* current = table[i];
        while (current) {
            HashNode* temp = current;
            current = current->next;
            free(temp->kmer);
            free(temp);
        }
        table[i] = NULL;
    }
}
