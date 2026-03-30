#ifndef SKIPLIST_H
#define SKIPLIST_H

typedef struct SkipNode {
    char* sequence;
    int score;
    struct SkipNode** forward;
} SkipNode;

void insert_skiplist(const char *seq, const char *species, int score);
void display_top_matches(int k);
void free_skiplist();
void init_skiplist();

#endif