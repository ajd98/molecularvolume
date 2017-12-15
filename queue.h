#ifndef QUEUE_H
#define QUEUE_H

// encapsulate grid indices
struct Triple
{
    int x, y, z;
};

// fifo data structure
struct Queue 
{
    int front, rear, size;
    unsigned int capacity;
    struct Triple* items;
};

struct Queue* newQueue(int capacity);
int isFull(struct Queue* queue);
int isEmpty(struct Queue* queue);
int append(struct Triple coords, struct Queue* queue);
struct Triple pop(struct Queue* queue);

#endif
