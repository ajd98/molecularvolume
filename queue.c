#include "queue.h"

// Make a new queue with maximum capacity of ``capacity`` Triples
int newQueue(int capacity)
{
    struct Queue* queue;
    queue->front = 0;
    queue->size = 0;
    queue->rear = 0;
    queue->capacity = capacity;
    queue->items = (struct Triple*)malloc(capacity*sizeof(struct Triple));
}

// Check if the queue is full
int isFull(struct Queue* queue)
{
    if (queue->capacity == queue->size){
        return 1;
    } else {
        return 0;
    }
}

// Add a Triple to the queue
int enqueue(struct Triple coords, struct Queue* queue)
{
    if isFull(queue){
        return 1; // error; queue is full
    }
    // Add the item to the rear of the queue
    // wrapping around to the beginning of the array if need be
    queue->items[queue->rear % queue->capacity]  = coords;
    queue->rear = (queue->rear+1) % queue->capacity;
    queue->size += 1;
    return 0;
}

int isEmpty(struct Queue* queue)
{
    if (queue->size == 0){
        return 1;
    } else {
        return 0;
    }
    
}

struct Triple dequeue(struct Queue* queue)
{
    if isEmpty(queue){
        return NULL; // queue is empty; return a null pointer
    }
    struct Triple item = queue->items[queue->front];
    queue->front = (queue->front+1)%queue->capacity;
    queue->size -= 1;
    return item;
}

int main()
{
    struct Queue queue = newQueue(10);
    struct Triple item1;
    item1.x = 1;
    item1.y = 1;
    item1.z = 1;
    enqueue(item1, queue);
}
