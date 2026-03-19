
#pragma once

#ifndef TURBINECORE_STACK_H
#define TURBINECORE_STACK_H

#include "wind_turbine_core/ProjectDefines.h"

namespace turbine_core {

    namespace topology {

        namespace stack {

            template< typename T >
            struct Stack {
                unsigned int capacity;
                int top;
                T ** node;
            };

            template< typename T >
            HOST_DEVICE_PREFIX Stack<T> * createStack(unsigned capacity){
                auto * stack = (Stack<T>*)malloc(sizeof(Stack<T>));
                stack->capacity = capacity;
                stack->top = -1;
                stack->node = (T**)malloc(capacity * sizeof(T*)); // NOLINT(bugprone-sizeof-expression)

                return stack;
            }

            template< typename T >
            HOST_DEVICE_PREFIX int isFull(Stack<T> * stack){
                return (stack->top == (stack->capacity - 1));
            }

            template< typename T >
            HOST_DEVICE_PREFIX int isEmpty(Stack<T> * stack){
                return (stack->top == -1);
            }

            template< typename T >
            HOST_DEVICE_PREFIX void push(Stack<T> * stack, T * item){
                if(isFull(stack)){
                    return;
                }
                stack->top += 1;
                stack->node[stack->top] = item;
            }

            template< typename T >
            HOST_DEVICE_PREFIX T * peek(Stack<T> * stack){
                if(isEmpty(stack)){
                    return nullptr;
                }

                return stack->node[stack->top];
            }

            template< typename T >
            HOST_DEVICE_PREFIX T * pop(Stack<T> * stack){
                if(isEmpty(stack)){
                    return nullptr;
                }
                T * node =  stack->node[stack->top];
                stack->top -= 1;

                return node;
            }

        }

    }

}

#endif //TURBINECORE_STACK_H
