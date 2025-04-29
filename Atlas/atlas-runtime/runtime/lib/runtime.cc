#include "runtime.h"
#include "bks_ctx.h"
#include "helpers.h"

atlas::BksContext *bks_ctx = nullptr;

void runtime_init() {
    if (!bks_ctx)
        bks_ctx = new atlas::BksContext();
}

void runtime_exit() {
    delete bks_ctx;
    bks_ctx = nullptr;
}

void *runtime_fetch(const void *object, int object_size) {
    BARRIER_ASSERT(bks_ctx != nullptr);
    return bks_ctx->Fetch(object, object_size);
}

bool runtime_read(void *dst, const void *object, int object_size) {
    BARRIER_ASSERT(bks_ctx != nullptr);
    return bks_ctx->Read(dst, object, object_size);
}

// static FORCE_INLINE void timer_start(unsigned *cycles_high_start,
//     unsigned *cycles_low_start) {
// asm volatile("xorl %%eax, %%eax\n\t"
// "CPUID\n\t"
// "RDTSC\n\t"
// "mov %%edx, %0\n\t"
// "mov %%eax, %1\n\t"
// : "=r"(*cycles_high_start), "=r"(*cycles_low_start)::"%rax",
// "%rbx", "%rcx", "%rdx");
// }

// static FORCE_INLINE void timer_end(unsigned *cycles_high_end,
//   unsigned *cycles_low_end) {
// asm volatile("RDTSCP\n\t"
// "mov %%edx, %0\n\t"
// "mov %%eax, %1\n\t"
// "xorl %%eax, %%eax\n\t"
// "CPUID\n\t"
// : "=r"(*cycles_high_end), "=r"(*cycles_low_end)::"%rax", "%rbx",
// "%rcx", "%rdx");
// }

// static FORCE_INLINE uint64_t get_elapsed_cycles(unsigned cycles_high_start,
//     unsigned cycles_low_start,
//     unsigned cycles_high_end,
//     unsigned cycles_low_end) {
// uint64_t start, end;
// start = ((static_cast<uint64_t>(cycles_high_start) << 32) | cycles_low_start);
// end = ((static_cast<uint64_t>(cycles_high_end) << 32) | cycles_low_end);
// return end - start;
// }