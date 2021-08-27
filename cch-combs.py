# Combinations: {A, A+1, ..., B} choose K
#
# Decomposes each combination into two halves (plus one element between them if K is odd).
# Each combination is a tuple (X, Y, ...), where:
# X is the first/lowest element, Y the last/highest;
# - for a single element, a combination is the 2-tuple (X, X)
# - for an even count, a 4-tuple (X, Y, C, D), where C and D are combinations
# - for an odd count, a 5-tuple (X, Y, C, M, D), where M is a choose-1 combination
# C, D and M follow the same pattern, recursively
#
# This awkward representation mimics the way combinations are stored in cch-detect.cpp,
# as indices of two combinations of half as many elements (plus one if odd); it's not
# useful as such. See combs2() for the more logical approach. See combs1() for what
# actually went into cch-detect.cpp

def combs(A, B, K):
	if K == 1:
		print(f'{A}..{B} choose one')
		return [(M, M) for M in range(A, B + 1)]

	J = K // 2
	if K % 2 == 0:
		print(f'{A}..{B} choose {K}, even, two halves of {J}')
		return [(C[0], D[1], C, D) for C in combs(A, B - J, J) for D in combs(C[1] + 1, B, J)]
	else:
		print(f'{A}..{B} choose {K}, odd, two halves of {J} plus one')
		return [(C[0], D[1], C, M, D) for C in combs(A, B - J - 1, J) for M in combs(C[1] + 1, B - J, 1) for D in combs(M[1] + 1, B, J)]


# Generate choose-K combinations that have x as the first element and y as the last one.
# Build them into a singly-linked list and memorize the list head for a given (K, x, y).
# Build larger combinations by joining two halves (plus a single element if odd).

def combs1(N, K):
	# make and index {0, ..., N-1} choose 1
	res = [('single', x) for x in range(N)]
	nxt = [None] * N
	res1 = []
	mem = {}  # K = 2 and 3 are special-cased below, so no need for choose-1 to be memorized

	def comb_at(x, y, K, final):
		#print(f'comb_at({x}, {y}, {K})')
		#input('...')

		if K == 1:
			assert(x == y)
			return x

		if (K, x, y) in mem:
			return mem[(K, x, y)]

		pos = None
		ret = None
		if K == 2:
			if final:
				res1.append((x, y))
			else:
				res.append((x, y))
				mem[(K, x, y)] = ret = pos = len(nxt)
				nxt.append(None)
		elif K == 3:
			for m in range(x + 1, y):
				if final:
					res1.append((x, m, y))
				else:
					res.append((x, m, y))
					if pos is None:
						mem[(K, x, y)] = ret = pos = len(nxt)
					else:
						nxt[pos] = len(nxt)
						pos = len(nxt)
					nxt.append(None)
		else:
			J = K // 2
			if K % 2 == 0:
				for y1 in range(x + J - 1, y - J + 1):
					i0 = comb_at(x, y1, J, False)
					for x1 in range(y1 + 1, y - J + 2):
						j0 = comb_at(x1, y, J, False)
						i = i0
						while i is not None:
							j = j0
							while j is not None:
								if final:
									res1.append((i, j))
								else:
									res.append((i, j))
									if pos is None:
										mem[(K, x, y)] = ret = pos = len(nxt)
									else:
										nxt[pos] = len(nxt)
										pos = len(nxt)
									nxt.append(None)
								j = nxt[j]
							i = nxt[i]
			else:
				for y1 in range(x + J - 1, y - J):
					i0 = comb_at(x, y1, J, False)
					for x1 in range(y1 + 2, y - J + 2):
						j0 = comb_at(x1, y, J, False)
						i = i0
						while i is not None:
							for m in range(y1 + 1, x1):
								j = j0
								while j is not None:
									if final:
										res1.append((i, m, j))
									else:
										res.append((i, m, j))
										if pos is None:
											mem[(K, x, y)] = ret = pos = len(nxt)
										else:
											nxt[pos] = len(nxt)
											pos = len(nxt)
										nxt.append(None)
									j = nxt[j]
							i = nxt[i]

		return ret

	for x in range(N - K + 1):
		for y in range(x + K - 1, N):
			comb_at(x, y, K, True)

	z = len(res)
	res.extend(res1)

	return res, z
	

def combs2(A, B, K):
	if K == 1:
		return [(X,) for X in range(A, B + 1)]

	J = K // 2
	if K % 2 == 0:
		return [C + D for C in combs2(A, B - J, J) for D in combs2(C[-1] + 1, B, J)]
	else:
		return [C + M + D for C in combs2(A, B - J - 1, J) for M in combs2(C[-1] + 1, B - J, 1) for D in combs2(M[-1] + 1, B, J)]


def extract_comb(X):
	if len(X) == 2:
		return (X[0],)
	elif len(X) == 4:
		A, B, C, D = X
		return extract_comb(C) + extract_comb(D)
	elif len(X) == 5:
		A, B, C, M, D = X
		return extract_comb(C) + extract_comb(M) + extract_comb(D)

def extract_comb1(i, R):
	if R[i][0] == 'single':
		return (R[i][1],)
	else:
		if len(R[i]) == 2:
			return extract_comb1(R[i][0], R) + extract_comb1(R[i][1], R)
		else:
			return extract_comb1(R[i][0], R) + extract_comb1(R[i][1], R) + extract_comb1(R[i][2], R)

def main():
	N, K = 21, 11
	A, B = 0, N - 1
	
	#L1 = combs(A, B, K)
	L2 = combs2(A, B, K)
	#L3 = [extract_comb(X) for X in L1]
	#for X in L1:
	#	print(X)
	#for X in L3:
	#	print(X)
	#print(L2 == L3)
	print(len(L2))

	R, z = combs1(N, K)
	print(len(R) - z)
	for (i, x) in enumerate(R):
		s = '*' if i >= z else ''
		print(f'{i}{s}: {x}')

	for i in range(z, len(R)):
		print(extract_comb1(i, R))

	print(len(R) - z)



if __name__ == "__main__":
	main()
