#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "suffix_tree.h"

int get_index(char c) {
    switch (c) {
        case 'A': return 0;
        case 'T': return 1;
        case 'C': return 2;
        case 'G': return 3;
        default: return -1;
    }
}

SuffixTreeNode* create_node() {
    SuffixTreeNode* node = (SuffixTreeNode*)malloc(sizeof(SuffixTreeNode));
    for (int i = 0; i < 4; i++) {
        node->children[i] = NULL;
    }
    node->start = -1;
    node->end = -1;
    return node;
}

SuffixTreeNode* root = NULL;

void insert_suffix(const char* text, int start_index) {
    SuffixTreeNode* current = root;

    for (int i = start_index; text[i] != '\0'; i++) {
        int idx = get_index(text[i]);
        if (idx == -1) continue;

        if (current->children[idx] == NULL) {
            current->children[idx] = create_node();
        }

        current = current->children[idx];
    }
}

void build_suffix_tree(const char* text) {
    if (root != NULL) {
        free_suffix_tree();
    }

    root = create_node();

    int n = strlen(text);
    for (int i = 0; i < n; i++) {
        insert_suffix(text, i);
    }
}

int search_pattern(const char* pattern) {
    if (root == NULL) {
        printf("Suffix Tree not built.\n");
        return 0;
    }

    SuffixTreeNode* current = root;

    for (int i = 0; pattern[i] != '\0'; i++) {
        int idx = get_index(pattern[i]);
        if (idx == -1) return 0;

        if (current->children[idx] == NULL) {
            return 0;
        }

        current = current->children[idx];
    }

    return 1; 
}

void free_node(SuffixTreeNode* node) {
    if (node == NULL) return;

    for (int i = 0; i < 4; i++) {
        free_node(node->children[i]);
    }

    free(node);
}

void free_suffix_tree() {
    free_node(root);
    root = NULL;
}