from filtros.chebyshev import chebyshevFilter


def main():
	order = 5
	fc = 10000
	rp = 3

	filtro = chebyshevFilter(order, fc, rp)
	filtro.calculate_filters_values()


if __name__ == '__main__':
	main()