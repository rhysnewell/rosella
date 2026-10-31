// The host's coverage constants. Metal has no doubles.
struct Terms {
    float weight, heaviest, epsilon, min_var_epsilon;
};

// `scale` is 1 / (deviation * sqrt 2), so a cdf costs a product and not a division.
struct Moment {
    float mean, variance, scale, log_variance;
};

// Metal has no erfc. Press et al. (2007) fit it to a relative error under 1.2e-7 everywhere.
static inline float erfc_of(float z) {
    float t = 1.0f / (1.0f + 0.5f * fabs(z));
    float poly = -1.26551223f + t * (1.00002368f + t * (0.37409196f + t * (0.09678418f
        + t * (-0.18628806f + t * (0.27886807f + t * (-1.13520398f + t * (1.48851587f
        + t * (-0.82215223f + t * 0.17087277f))))))));
    float tail = t * exp(-z * z + poly);
    return z >= 0 ? tail : 2.0f - tail;
}

static inline float normal_cdf(Moment m, float x) {
    return fmin(0.5f * erfc_of((m.mean - x) * m.scale), 1.0f);
}

// The host measures one pair from either end, and both have to round alike.
static float overlap(Moment a, Moment b, device const Terms* terms) {
    if (b.variance < a.variance || (b.variance == a.variance && b.mean < a.mean)) {
        Moment swap = a;
        a = b;
        b = swap;
    }
    float k1, k2;
    if (fabs(a.variance - b.variance) < terms->min_var_epsilon) {
        k1 = k2 = (a.mean + b.mean) / 2;
    } else {
        float gap = (a.mean - b.mean) * (a.mean - b.mean)
            - (a.variance - b.variance) * (b.log_variance - a.log_variance);
        float tmp = sqrt(a.variance * b.variance * gap);
        float across = 1 / (a.variance - b.variance);
        k1 = (tmp - a.mean * b.variance + b.mean * a.variance) * across;
        k2 = -(tmp + a.mean * b.variance - b.mean * a.variance) * across;
    }
    if (k1 > k2) {
        float swap = k1;
        k1 = k2;
        k2 = swap;
    }
    Moment narrow = a.variance > b.variance ? b : a;
    Moment wide = a.variance > b.variance ? a : b;
    if (k1 == k2) return fabs(normal_cdf(narrow, k1) - normal_cdf(wide, k1));
    return fabs(normal_cdf(narrow, k2) - normal_cdf(narrow, k1)
                + normal_cdf(wide, k1) - normal_cdf(wide, k2));
}

static inline float composition_only(float dot, float mine, float theirs) {
    float total = mine + theirs;
    float distance = 0.0f;
    if (total != 0.0f) distance = (total - 2.0f * dot) / total;
    return isnan(distance) ? 2.0f : distance;
}

static inline bool absent(float x_mean, float y_mean, float mine_bar, float their_bar, float epsilon) {
    return x_mean - epsilon <= mine_bar && y_mean - epsilon <= their_bar;
}

// Each overlap only adds to the sum, so once the sum alone puts the pair past `worst` the rest
// cannot bring it back. The composition floor is the same bound before any sample is read.
static float mixed(float dot, float mine, float theirs, uint a, uint b, float worst,
                   device const float* moments, device const float* presence, int stride, int samples,
                   device const Terms* terms) {
    float total = mine + theirs;
    float composition = total == 0 ? 0.0f : (total - 2 * dot) / total;
    if (isnan(composition)) composition = 2;
    if (composition > 0 && (1 - terms->heaviest) * composition > worst) return INFINITY;
    float epsilon = terms->epsilon, mine_bar = presence[a], their_bar = presence[b];
    int scored = 0;
    for (int s = 0; s < samples; ++s) {
        device const float* sample = moments + (long)s * 4 * stride;
        scored += !absent(sample[a], sample[b], mine_bar, their_bar, epsilon);
    }
    float weight = terms->weight < 0 ? (float)scored / ((float)scored + 1) : terms->weight;
    float sum = 0, enough = INFINITY;
    if (scored > 0 && weight > 0 && isfinite(worst)) {
        enough = (worst - (1 - weight) * composition) / weight * scored * (1 + 1e-5f) + 1e-5f;
    }
    for (int s = 0; s < samples && scored > 0; ++s) {
        device const float* sample = moments + (long)s * 4 * stride;
        float x_mean = sample[a], y_mean = sample[b];
        if (absent(x_mean, y_mean, mine_bar, their_bar, epsilon)) continue;
        Moment x = {x_mean, sample[stride + a], sample[2 * stride + a], sample[3 * stride + a]};
        Moment y = {y_mean, sample[stride + b], sample[2 * stride + b], sample[3 * stride + b]};
        float shared = overlap(x, y, terms);
        sum += shared < epsilon ? epsilon : (shared > 1 - epsilon ? 1 - epsilon : shared);
        if (sum > enough) return INFINITY;
    }
    float coverage = epsilon;
    if (scored > 0) {
        coverage = sum / scored;
        if (isnan(coverage)) coverage = 1;
    }
    float distance = weight * coverage + (1 - weight) * composition;
    return isnan(distance) ? 1.0f : distance;
}

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
                    device float* distances [[buffer(8)]], device uint* indices [[buffer(9)]],
                    device const uint* query_rows [[buffer(10)]], device const float* moments [[buffer(11)]],
                    device const float* presence [[buffer(12)]], constant int& stride [[buffer(13)]],
                    constant int& samples [[buffer(14)]], device const float* terms [[buffer(15)]], INDEX) {
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
    uint query = samples > 0 ? query_rows[row] : 0u;
    device const float* line = dots + (long)row * cols;
    for (int start = 0; start < cols; start += 32) {
        int col = start + lane;
        uint at = (uint)(first + col);
        bool live = col < cols && at != own;
        float distance = INFINITY;
        if (live) {
            distance = samples > 0
                ? mixed(line[col], mine, base_variance[at], query, at, worst, moments, presence, stride,
                        samples, (device const Terms*)terms)
                : composition_only(line[col], mine, base_variance[at]);
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
