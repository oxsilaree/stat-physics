#include "parameters.h"
#include "functions.h"

int* getNeighbor(int i, int j, int pos)
{ // 0 -> left, 1 -> right, 2 -> up, 3 -> down, 4 -> up(NNN), 5 -> down(NNN)
    static int r[2] = {0,0};
    if (pos == 0) { // left
        r[1] = j;
        if (i == 0){
            r[0] = i + (L-1);
        } else {
            r[0] = i - 1;  
        }
    }
    if (pos == 1) { // right
        r[1] = j;
        if (i == L - 1){
            r[0] = i - (L-1);
        } else {
            r[0] = i + 1;  
        }
    }
    if (pos == 2) { // up
        r[0] = i;
        if (j == 0){
            r[1] = j + (L-1);
        } else {
            r[1] = j - 1;  
        }
    }
    if (pos == 3) { // down
        r[0] = i;
        if (j == L - 1){
            r[1] = j - (L-1);
        } else {
            r[1] = j + 1;  
        }
    }
    if (pos == 4) { // up NNN
        r[0] = i;
        if (j <= 1){
            r[1] = j + (L-2);
        } else {
            r[1] = j - 2;  
        }
    }
    if (pos == 5) { // down NNN
        r[0] = i;
        if (L-2 <= j){
            r[1] = j - (L-2);
        } else {
            r[1] = j + 2;  
        }
    }
    return r;
}
