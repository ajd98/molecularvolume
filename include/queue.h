#ifndef QUEUE_H
#define QUEUE_H

// encapsulate grid indices
struct Triple
{
    unsigned short int x, y, z;
};

// fifo data structure
struct Queue 
{
    unsigned long front, rear, size;
    unsigned long capacity;
    struct Triple* items;
};

struct Queue* newQueue(unsigned long capacity);
void delQueue(struct Queue* queue);
int isFull(struct Queue* queue);
int isEmpty(struct Queue* queue);
int append(struct Triple coords, struct Queue* queue);
struct Triple pop(struct Queue* queue);

#endif
