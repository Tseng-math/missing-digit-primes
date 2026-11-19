// Directed fixed-point verification of a shifted pure-prime word norm.
// stdin supplies an upper table for the one-digit mask on a 2^17 grid
// in units of 2^36, produced by pure.py.
#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

using u64 = std::uint64_t;
using wide = unsigned __int128;
constexpr u64 scale = u64(1) << 36;
constexpr u64 mesh = u64(1) << 17;
constexpr int max_depth = 50;
std::vector<u64> mask(mesh);
std::vector<u64> modulus, multiplier;
std::vector<int> freeze_depth;
u64 prime, cap_num, cap_den, leaves = 0;
int max_leaf_depth = 0;

u64 up(wide a, wide b) {
    return (a + b - 1) / b;
}

// An array entry is an upper bound for the accumulated weights of the
// corresponding residue in the current pure-prime grid. The midpoint
// shift is (2*prefix+1)/2^(depth+1). Contract the highest p-adic digit.
std::vector<u64> fold(const std::vector<u64>& buf, int position,
                      u64 prefix, int depth) {
    const u64 q = modulus[position];
    const u64 c = multiplier[position];
    const wide denominator = wide(q) << (depth + 1);
    const u64 guard = wide(19) * c >= denominator ? scale :
                      up(wide(19) * scale, 2 * mesh)
                      + up(wide(19) * scale * c, denominator);
    const u64 step = (c % q) * mesh;
    const u64 period = q * mesh;
    const u64 length = q / prime;
    const u64 rounding_offset =
        ((wide(c) * (2 * prefix + 1) * mesh + denominator / 2) >> (depth + 1))
        % period;
    const u64* weights = buf.empty() ? nullptr : buf.data();
    std::vector<u64> out(length, 0);
    if (guard >= scale) {
        for (u64 i = 0; i < prime; ++i)
            for (u64 j = 0; j < length; ++j)
                out[j] += weights ? weights[i * length + j] : scale;
        return out;
    }
    for (u64 i = 0; i < prime; ++i) {
        u64 phase = (rounding_offset + wide(i * length) * step) % period;
        for (u64 j = 0; j < length; ++j) {
            const u64 idx = (phase / q) & (mesh - 1);
            const u64 upper = std::min(scale, mask[idx] + guard);
            const wide product = wide(upper) * (weights ? weights[i * length + j] : scale);
            out[j] += (product + scale - 1) >> 36;
            phase += step;
            if (phase >= period) phase -= period;
        }
    }
    return out;
}

void certify(const std::vector<u64>& frozen, int pos, u64 prefix, int depth) {
    // At the first depth where a low-factor interval guard is at most
    // 1/4096, reuse its directed upper bound throughout all descendants.
    // These factors can be peeled in order since their slopes increase.
    const std::vector<u64>* initial = &frozen;
    std::vector<u64> new_frozen;
    const int end = modulus.size();
    while (pos < end && depth >= freeze_depth[pos]) {
        new_frozen = fold(*initial, pos, prefix, depth);
        initial = &new_frozen;
        ++pos;
    }
    const std::vector<u64>* next = initial;
    std::vector<u64> active;
    for (int j = pos; j < end; ++j) {
        active = fold(*next, j, prefix, depth);
        next = &active;
    }
    if (wide(next->at(0)) * cap_den <= wide(scale) * cap_num) {
        ++leaves;
        max_leaf_depth = std::max(max_leaf_depth, depth);
        return;
    }
    if (depth == max_depth) {
        std::cerr << "uncertified interval prefix " << prefix << ", depth "
                  << depth << ", upper integer " << next->at(0) << '\n';
        std::exit(1);
    }
    certify(*initial, pos, 2 * prefix, depth + 1);
    certify(*initial, pos, 2 * prefix + 1, depth + 1);
}

int main() {
    u64 power;
    if (!(std::cin >> prime >> power >> cap_num >> cap_den) ||
        (prime != 2 && prime != 5) || power == 0 || power > 24 ||
        cap_den == 0) return 1;
    for (u64& node : mask) {
        if (!(std::cin >> node) || node > scale) return 1;
    }
    u64 q = 1, c = 1;
    for (u64 i = 0; i < power; ++i) q *= prime;
    if (q > (u64(1) << 24)) return 1;
    for (u64 i = 0; i < power; ++i, q /= prime, c *= 10 / prime) {
        modulus.push_back(q);
        multiplier.push_back(c);
        int depth = 0;
        while (wide(19) * c * 4096 > (wide(q) << (depth + 1))) ++depth;
        freeze_depth.push_back(depth);
    }
    certify({}, 0, 0, 0);
    std::cout << prime << ' ' << power << ": norm <= " << cap_num << '/' << cap_den
              << " on " << leaves << " intervals, max depth " << max_leaf_depth << '\n';
}
