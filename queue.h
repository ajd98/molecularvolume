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

int newQueue(int capacity);
int isFull(struct Queue* queue);
int enqueue(struct Triple coords, struct Queue* queue);
struct Triple dequeue(struct Queue* queue);

#endif
