from filtros.chebyshev import chebyshevFilter


def main():
	order = 5
	fc = 10000
	rp = 3

	filtro = chebyshevFilter(order, fc, rp)
	print(filtro.response)

if __name__ == '__main__':
	main()
