import numpy as np

def edit_nested_plain(S1, P1, S2):
    n = len(S1)
    m = len(S2)

    # Initialize the DP table
    DP = np.full((n+2, n+2, m+2, m+2), np.inf)

    # Scoring scheme
    wm = 1  # base-mismatch
    wd = 1  # base-deletion
    wa = 1  # arc-altering
    wb = 1  # arc-breaking
    wr = 1  # arc-removing

    def mismatch(i, j):
        return 1 if S1[i - 1] != S2[j - 1] else 0

    # Base case initialization
    for i in range(1, n+1):
        for j in range(1, m+1):
            DP[i, i - 1, j, j] = (j - j + 1) * wd

    # Main loop
    for k in range(1, n+1):
        for i in range(1, n - k + 2):
            i_prime = i + k - 1
            for j in range(1, m+1):
                for j_prime in range(j, m+1):

                    if (i, i_prime) in P1:
                        DP[i, i_prime, j, j_prime] = min(
                            DP[i+1, i_prime-1, j+1, j_prime-1] + wb + (mismatch(i, j) + mismatch(i_prime, j_prime)) * wm if j < j_prime else np.inf,
                            DP[i+1, i_prime-1, j, j_prime-1] + wa + mismatch(i_prime, j_prime) * wm,
                            DP[i+1, i_prime-1, j+1, j_prime] + wa + mismatch(i, j) * wm,
                            DP[i+1, i_prime-1, j, j_prime] + wr,
                            DP[i, i_prime, j, j_prime-1] + wd,
                            DP[i, i_prime, j+1, j_prime] + wd,
                        )
                    else:
                        if i_prime not in [u[1] for u in P1]:
                            DP[i, i_prime, j, j_prime] = min(
                                DP[i, i_prime-1, j, j_prime-1] + mismatch(i_prime, j_prime) * wm,
                                DP[i, i_prime-1, j, j_prime] + wd,
                                DP[i, i_prime, j, j_prime-1] + wd,
                            )
                        else:
                            min_val = np.inf
                            for j_double_prime in range(j, j_prime+1):
                                cur_val = DP[i, [u[0] for u in P1 if u[1] == i_prime][0]-1, j, j_double_prime-1] + \
                                          DP[[u[0] for u in P1 if u[1] == i_prime][0], i_prime, j_double_prime, j_prime]
                                min_val = min(min_val, cur_val) 
                                DP[i, i_prime, j, j_prime] = min_val
    return DP[1, n, 1, m]                                
S1 = "GCAUC"
P1 = [(1, 4), (2, 3)] # Nested structure represented as pairs of indices
S2 = "GGCAU"
edit_distance = edit_nested_plain(S1, P1, S2)
print("Edit distance:", edit_distance)