#define NEAREST_SLOTS 4

static inline bool ahead(float d, uint at, float e, uint bt) {
    return d < e || (d == e && at < bt);
}

static inline uint ballot(bool flag) {
    return (uint)(simd_vote::vote_t)simd_ballot(flag);
}

// A SIMD group per query holds its k nearest in registers, slot s in lane s % 32, so an insert is
// a few shuffles and not a walk through memory. Ties go to the lower index as the host lists do.
kernel void nearest(device const float* dots [[buffer(0)]], constant int& rows [[buffer(1)]],
                    constant int& cols [[buffer(2)]], constant int& first [[buffer(3)]],
                    device const uint* selves [[buffer(4)]], device const float* query_variance [[buffer(5)]],
                    device const float* base_variance [[buffer(6)]], constant int& k [[buffer(7)]],
                    device float* distances [[buffer(8)]], device uint* indices [[buffer(9)]], INDEX) {
    int row = FLAT / 32;
    uint lane = local.x & 31;
    if (row >= rows) return;
    device float* best = distances + (long)row * k;
    device uint* held = indices + (long)row * k;
    float near[NEAREST_SLOTS];
    uint near_at[NEAREST_SLOTS];
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        int slot = r * 32 + lane;
        near[r] = slot < k ? best[slot] : INFINITY;
        near_at[r] = slot < k ? held[slot] : 0xffffffffu;
    }
    int last_reg = (k - 1) / 32;
    ushort last_lane = (k - 1) % 32;
    float worst = 0.0f;
    uint worst_at = 0;
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        float d = simd_shuffle(near[r], last_lane);
        uint a = simd_shuffle(near_at[r], last_lane);
        if (r == last_reg) { worst = d; worst_at = a; }
    }
    float mine = query_variance[row];
    uint own = selves[row];
    device const float* line = dots + (long)row * cols;
    for (int start = 0; start < cols; start += 32) {
        int col = start + lane;
        uint at = (uint)(first + col);
        bool live = col < cols && at != own;
        float distance = INFINITY;
        if (live) {
            float total = mine + base_variance[at];
            distance = 0.0f;
            if (total != 0.0f) distance = (total - 2.0f * line[col]) / total;
            if (isnan(distance)) distance = 2.0f;
        }
        uint wanted = ballot(live && ahead(distance, at, worst, worst_at));
        while (wanted) {
            ushort from = ctz(wanted);
            wanted &= wanted - 1;
            float d = simd_shuffle(distance, from);
            uint a = simd_shuffle(at, from);
            if (!ahead(d, a, worst, worst_at)) continue;
            int place = 0;
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                bool before = r * 32 + (int)lane < k && ahead(near[r], near_at[r], d, a);
                place += popcount(ballot(before));
            }
            float up[NEAREST_SLOTS];
            uint up_at[NEAREST_SLOTS];
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                up[r] = simd_shuffle_up(near[r], 1);
                up_at[r] = simd_shuffle_up(near_at[r], 1);
                float wrapped = simd_shuffle(r > 0 ? near[r - 1] : 0.0f, 31);
                uint wrapped_at = simd_shuffle(r > 0 ? near_at[r - 1] : 0u, 31);
                if (lane == 0) { up[r] = wrapped; up_at[r] = wrapped_at; }
            }
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                int slot = r * 32 + lane;
                if (slot > place) { near[r] = up[r]; near_at[r] = up_at[r]; }
                else if (slot == place) { near[r] = d; near_at[r] = a; }
            }
            for (int r = 0; r < NEAREST_SLOTS; ++r) {
                float e = simd_shuffle(near[r], last_lane);
                uint b = simd_shuffle(near_at[r], last_lane);
                if (r == last_reg) { worst = e; worst_at = b; }
            }
        }
    }
    for (int r = 0; r < NEAREST_SLOTS; ++r) {
        int slot = r * 32 + lane;
        if (slot < k) { best[slot] = near[r]; held[slot] = near_at[r]; }
    }
}
