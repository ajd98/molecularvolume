#include "queue.h"
#include <stdio.h>
#include <stdlib.h>

// Make a new queue with maximum capacity of ``capacity`` Triples
struct Queue* newQueue(int capacity)
{
    struct Queue* queue = (struct Queue*)malloc(sizeof(struct Queue));
    queue->front = 0;
    queue->size = 0;
    queue->rear = 0;
    queue->capacity = capacity;
    queue->items = (struct Triple*)malloc(capacity*sizeof(struct Triple));
    return queue;
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

int isEmpty(struct Queue* queue)
{
    if (queue->size == 0){
        return 1;
    } else {
        return 0;
    }
    
}

// Add a Triple to the queue
int append(struct Triple coords, struct Queue* queue)
{
    if (isFull(queue)){
        return 1; // error; queue is full
    }
    // Add the item to the rear of the queue
    // wrapping around to the beginning of the array if need be
    printf("adding item to position %d\n", queue->rear % queue->capacity);
    queue->items[queue->rear % queue->capacity]  = coords;
    queue->rear = (queue->rear+1) % queue->capacity;
    queue->size += 1;
    return 0;
}

// Always check that the queue is not empty before calling this!
struct Triple pop(struct Queue* queue)
{
    printf("taking item from position %d\n", queue->front);
    struct Triple item = queue->items[queue->front];
    queue->front = (queue->front+1)%(queue->capacity);
    if (queue->size > 0){
        queue->size -= 1;
    }
    return item;
}

int main()
{
    struct Queue* queue = newQueue(3);
    struct Triple item0;
    item0.x = 0;
    item0.y = 0;
    item0.z = 0;

    struct Triple item1;
    item1.x = 1;
    item1.y = 1;
    item1.z = 1;

    struct Triple item2;
    item2.x = 2;
    item2.y = 2;
    item2.z = 2;

    struct Triple item3;
    item3.x = 3;
    item3.y = 3;
    item3.z = 3;

    printf("appending item0 to queue\n");
    append(item0, queue);
    printf("appending item1 to queue\n");
    append(item1, queue);
    printf("appending item2 to queue\n");
    append(item2, queue);

    printf("finished appending; now popping\n");
    struct Triple returneditem = pop(queue);
    printf("popped item.x, item.y, item.z: %d %d %d\n", returneditem.x, returneditem.y, returneditem.z);

    printf("appending item3 to queue\n");
    append(item3, queue);

    while (!isEmpty(queue)){
        printf("size of queue: %d", queue->size);
        returneditem = pop(queue);
        printf("popped item.x, item.y, item.z: %d %d %d\n", returneditem.x, returneditem.y, returneditem.z);
    }

    return 0;
}
