from filtros.chebyshev import chebyshevFilter
from filtros.butterworth import butterFilter
from filtros.bessel import besselFilter


def main():
	order = 3
	fc = 1000
	rp = 0.5

	filtro = butterFilter(order, fc)
	filtro.print_results()

if __name__ == '__main__':
	main()
