#ifndef GROUP_H
#define GROUP_H

#define PMAX 1080
#define CMAX 17

typedef int group_t;

extern unsigned P;
extern unsigned C;
extern double ReChar[CMAX][CMAX];
extern double ImChar[CMAX][CMAX];
extern double ReTr[PMAX];
extern double ImTr[PMAX];
extern group_t mult[PMAX][PMAX];
extern group_t inv[PMAX];
extern group_t conclass[PMAX];
extern group_t id;

void load_group(const char *);

#endif /* ndef GROUP_H */
