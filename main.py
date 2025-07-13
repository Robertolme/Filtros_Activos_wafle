from filtros.chebyshev import chebyshevFilter


def main():
	order = 5
	fc = 10000
	rp = 3

	filtro = chebyshevFilter(order, fc, rp)
	filtro.print_results()

if __name__ == '__main__':
	main()
