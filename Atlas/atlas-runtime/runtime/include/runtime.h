#pragma once

#ifdef __cplusplus
extern "C" {
#endif

void runtime_init();
void runtime_exit();

/* low-level apis */
void *runtime_fetch(const void *object, int object_size);
bool runtime_read(void *dst, const void *object, int object_size);

#ifdef __cplusplus
}
#endif


// typedef unsigned long uint64_t;
// static void timer_start(unsigned *cycles_high_start,
//     unsigned *cycles_low_start);
// static void timer_end(unsigned *cycles_high_end, unsigned *cycles_low_end);
// static uint64_t get_elapsed_cycles(unsigned cycles_high_start,
//     unsigned cycles_low_start,
//     unsigned cycles_high_end,
//     unsigned cycles_low_end);