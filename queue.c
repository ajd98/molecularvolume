#include "queue.h"

struct Triple
{
    int x, y, z;
};

struct Queue 
{
    int front, rear, size;
    unsigned int capacity;
    Triple* items;
}

int newQueue(int capacity)
{
    Queue* queue;
    queue->front = 0;
    queue->size = 0;
    queue->rear = 0;
    queue->capacity = capacity;
    queue->items = (Triple*)malloc(capacity*sizeof(Triple));
}

int isFull(Queue* queue)
{
    if (queue->capacity == queue->size){
        return 1;
    } else {
        return 0;
    }
}

int enqueue(Triple* coords, Queue* queue)
{
    if isFull(queue){
        return 1; // error; queue is full
    }
    // Add the item to the rear of the queue
    queue->items[queue->rear % queue->capacity]  = coords;
    queue->rear = (queue->rear+1) % queue->capacity;

    return 0;
}
