#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "trie.h"

static int char_to_index(char c) {
    switch (toupper((unsigned char)c)) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

TrieNode* create_trie() {
    TrieNode* node = (TrieNode*)malloc(sizeof(TrieNode));
    if (!node) return NULL;
    for (int i = 0; i < 4; ++i) {
        node->children[i] = NULL;
    }
    node->isEnd = 0;
    return node;
}

void insert_sequence(TrieNode* root, const char* seq) {
    if (!root || !seq) return;

    TrieNode* current = root;
    for (const char* p = seq; *p; ++p) {
        int idx = char_to_index(*p);
        if (idx < 0) continue; // ignore invalid chars
        if (current->children[idx] == NULL) {
            TrieNode* child = (TrieNode*)malloc(sizeof(TrieNode));
            if (!child) return; // allocation failure, abort (leak not handled in minimal C)
            for (int j = 0; j < 4; ++j) child->children[j] = NULL;
            child->isEnd = 0;
            current->children[idx] = child;
        }
        current = current->children[idx];
    }

    current->isEnd = 1;
}

int search_sequence(TrieNode* root, const char* seq) {
    if (!root || !seq) return 0;

    TrieNode* current = root;
    for (const char* p = seq; *p; ++p) {
        int idx = char_to_index(*p);
        if (idx < 0) return 0; // invalid letters in query, not found
        if (!current->children[idx]) return 0;
        current = current->children[idx];
    }

    return current->isEnd ? 1 : 0;
}

static void free_trie_node(TrieNode* node) {
    if (!node) return;
    for (int i = 0; i < 4; ++i) {
        if (node->children[i]) {
            free_trie_node(node->children[i]);
        }
    }
    free(node);
}

void free_trie(TrieNode* root) {
    free_trie_node(root);
}
