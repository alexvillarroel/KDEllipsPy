"""Prueba del módulo Fortran compilado con f2py."""
import numpy as np
import demo_f2py   # <- el .so que generó f2py

print("=== Módulo:", demo_f2py.__name__, "===")
print("Subrutinas disponibles:", [s for s in dir(demo_f2py) if not s.startswith("_")])

# --- 1) vecadd: suma de vectores ---
a = np.array([1.0, 2.0, 3.0])
b = np.array([10.0, 20.0, 30.0])
c = demo_f2py.vecadd(a, b)          # 'n' y 'c' los maneja f2py solo
print("\n[vecadd] a + b =", c)
assert np.allclose(c, [11.0, 22.0, 33.0])

# --- f2py autogenera la "firma" Python de cada subrutina: ---
print("\n[docstring autogenerado de vecadd]")
print(demo_f2py.vecadd.__doc__)

# --- 2) ellipse_patch: mini-mkstress ---
nx, ny = 40, 40
field = demo_f2py.ellipse_patch(nx, ny, 12.0, 6.0, 5.0, -1.0)
print("[ellipse_patch] shape:", field.shape,
      "| dentro:", field.max(), "| fuera:", field.min(),
      "| nº puntos dentro:", int((field > 0).sum()))

# Dibujo ASCII rápido de la elipse (sin depender de matplotlib)
print("\nMapa de la elipse (X = dentro):")
for j in range(ny - 1, -1, -2):
    row = "".join("X" if field[i, j] > 0 else "." for i in range(0, nx, 1))
    print(" ", row)

print("\nOK: Fortran -> .so -> NumPy funcionando.")
