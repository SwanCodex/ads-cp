#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "skiplist.h"

#define MAX 5   

typedef struct Node {
    char *sequence;
    int score; //percentage score
    struct Node *forward[MAX];
} Node;

static struct Node *header = NULL;

//create node 
struct Node* createNode(const char *seq, int score)
{
    struct Node *p = (struct Node *)malloc(sizeof(struct Node));

    p->sequence = (char *)malloc(strlen(seq) + 1);
    strcpy(p->sequence, seq);

    p->score = score;

    for (int i = 0; i < MAX; i++)
    {
        p->forward[i] = NULL;
    }

    return p;
}

int randlevel()
{
    int level = 1;

    while ((rand() % 2) && level < MAX)
    {
        level++;
    }

    return level;
}

//initialize skip list 
void init_skiplist()
{
    header = createNode("", -1);
    srand(time(NULL));
}

//insert function (DESCENDING ORDER based on score) 
void insert_skiplist(const char *seq, int score)
{
    if (header == NULL)
        init_skiplist();

    struct Node *update[MAX];
    struct Node *x = header;

    //traverse from top level to bottom 
    for (int i = MAX - 1; i >= 0; i--)
    {
        while (x->forward[i] != NULL &&
               x->forward[i]->score > score)   
        {
            x = x->forward[i];
        }
        update[i] = x;
    }

    x = x->forward[0];

    struct Node *nn = createNode(seq, score);
    int level = randlevel();

    for (int i = 0; i < level; i++)
    {
        nn->forward[i] = update[i]->forward[i];
        update[i]->forward[i] = nn;
    }
}

//display top K matches 
void display_top_matches(int k)
{
    if (header == NULL)
    {
        printf("Skip list is empty.\n");
        return;
    }

    struct Node *x = header->forward[0];
    int count = 0;

    printf("\nTop %d Matches:\n", k);

    while (x != NULL && count < k)
    {
        printf("%d. Score: %d%% | Sequence: %.50s...\n",
               count + 1,
               x->score,
               x->sequence);

        x = x->forward[0];
        count++;
    }

    if (count == 0)
    {
        printf("No matches found.\n");
    }
}

//display full skip list 
void display_skiplist()
{
    if (header == NULL)
    {
        printf("Skip list is empty.\n");
        return;
    }

    for (int i = MAX - 1; i >= 0; i--)
    {
        struct Node *x = header->forward[i];
        printf("Level %d: ", i);

        while (x != NULL)
        {
            printf("[%d%%] ", x->score);
            x = x->forward[i];
        }

        printf("\n");
    }
}

//free memory
void free_skiplist()
{
    if (header == NULL)
        return;

    struct Node *current = header->forward[0];

    while (current != NULL)
    {
        struct Node *temp = current;
        current = current->forward[0];

        free(temp->sequence);
        free(temp);
    }

    free(header);
    header = NULL;
}