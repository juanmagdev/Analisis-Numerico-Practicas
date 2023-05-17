from funciones import *

# Ejercicio 1
localizar_frontera_RK(lambda z: 1, 1)
localizar_frontera_RK(lambda z: 1 + z, 2)
localizar_frontera_RK(lambda z: 1 + z + z**2/2, 3)
localizar_frontera_RK(lambda z: 1 + z + z**2/2 + z*3/6, 4)
localizar_frontera_RK(lambda z: 1 + z + z**2/2 + z**3/6 + z**4/24, 5)

title("Fronteras RK")
legend(["RK1", "RK2", "RK3", "RK4", "RK5"])
show()