from math import floor, gcd, sqrt
from random import randrange
from numpy import mod

class EllipticCurve:
    def __init__(self, a: int, b: int, p: int):
        """Parámetros a, b sobre grupo p"""
        self.a = a
        self.b = b
        self.p = p
        
        # Si 4a^3 + 27b^2 ≠ 0 (mod p) 
        disc = (4 * (a ** 3) + 27 * (b ** 2)) % p
        if disc == 0:
            raise ValueError("La curva es singular")

class Point:
    def __init__(self, curve: EllipticCurve, x: int = None, y: int = None):
        """Inicializar un punto con sus coordenadas y al curva a la que pertenece. None representa el infinito"""
        self.curve = curve
        self.x = x
        self.y = y  

        if x is None and y is None:
            return
            
        # Comprobar si el punto esta en la curva
        if not self._is_on_curve():
            raise ValueError("El punto no está en la curva")

    
    def _is_on_curve(self) -> bool:
        """Comprobar si el punto esta en la curva"""
        if self.x is None or self.y is None:  # El punto en el infinito siempre pertenece
            return True
            
        #modulo p
        x = self.x % self.curve.p
        y = self.y % self.curve.p
        
        # Comprobar si cumple: y^2 = x^3 + ax + b (mod p)
        left = (y * y) % self.curve.p
        right = (pow(x, 3, self.curve.p) + 
                (self.curve.a * x) % self.curve.p + 
                self.curve.b) % self.curve.p
        return left == right
    
    def __eq__(self, other: 'Point') -> bool:
        """Comprobar si dos puntos son iguales"""
        if not isinstance(other, Point):
            return False
        if self.x is None or other.x is None:
            return self.x is None and other.x is None
        return (self.curve.a == other.curve.a and 
                self.curve.b == other.curve.b and 
                self.curve.p == other.curve.p and 
                self.x % self.curve.p == other.x % self.curve.p and 
                self.y % self.curve.p == other.y % self.curve.p)
    
    def __neg__(self) -> 'Point':
        """Devolver el punto negado"""
        if self.x is None:  # El punto en el infinito es su propio inverso
            return self
        return Point(self.curve, self.x, (-self.y) % self.curve.p)
    
    def den_inv(self, other: 'Point') -> int:
        if self.x is None or other.x is None:
            raise TypeError("No se puede calcular de un punto en el infinito")
        if self == other:
            den = (2 * self.y) % self.curve.p
        else:
            den = (other.x - self.x) % self.curve.p

        return pow(den, -1, self.curve.p)
    
    def __add__(self, other: 'Point') -> 'Point':
        if not isinstance(other, Point) or self.curve.p != other.curve.p:
            raise TypeError("Los puntos deben pertenecer a la misma curva")
        
        if self.x is None: return other
        if other.x is None: return self
        if self == -other: return Point(self.curve)
            
        try:
            if self == other:
                if self.y == 0:
                    return Point(self.curve)
                num = (3 * pow(self.x, 2, self.curve.p) + self.curve.a) % self.curve.p
                den = (2 * self.y) % self.curve.p
                if den == 0:
                    raise ZeroDivisionError(den)
                try:
                    slope = (num * pow(den, -1, self.curve.p)) % self.curve.p
                except ValueError:
                    raise ZeroDivisionError(den)
            else:
                num = (other.y - self.y) % self.curve.p
                den = (other.x - self.x) % self.curve.p
                if den == 0:
                    raise ZeroDivisionError(den)
                try:
                    slope = (num * pow(den, -1, self.curve.p)) % self.curve.p
                except ValueError:
                    raise ZeroDivisionError(den)
                
            x3 = (pow(slope, 2, self.curve.p) - self.x - other.x) % self.curve.p
            y3 = (-self.y + slope * (self.x - x3)) % self.curve.p
            
            return Point(self.curve, x3, y3)
            
        except ZeroDivisionError as e:
            raise ZeroDivisionError(e.args[0])
    
    def __sub__(self, other: 'Point') -> 'Point':
        """Resta de puntos"""
        return self + (-other)

    def __mul__(self, scalar: int) -> 'Point':
        """Multiplicacion punto por escalar usando el double-and-add"""
        if not isinstance(scalar, int):
            raise TypeError("El escalar tiene que ser un entero")
            
        # Casos especiales
        if scalar < 0:
            return (-self) * (-scalar) # El punto negado por escalar en positivo
        if scalar == 0:
            return Point(self.curve)  # Multiplicar por 0 devuelve el punto en el infinito
        if self.x is None:
            return self  # Multiplicar el punto en el infinito da él mismo
            
        # Double-and-add algorithm
        result = Point(self.curve)  # Empezamos con el punto en el infinito
        current = self
        while scalar:
            if scalar & 1:  # Si el bit es 1, se suma al resultado
                result += current
            current += current  # Se duplica el punto
            scalar >>= 1  # Siguiente bit
            
        return result
    
    def __rmul__(self, scalar: int) -> 'Point':
        """Multiplicación por la derecha"""
        return self * scalar

    def __str__(self) -> str:
        """Punto a string"""
        if self.x is None:
            return "Point(∞)"
        return f"Point({self.x}, {self.y})"
    
# Método pollardRho de descomposición de números en factores
def pollardRho(n):
    A = B = randrange(start=2, stop= n - 1, step=1)
    while True:
        A = (A*A + 1) % n
        B = (B*B + 1) % n
        B = (B*B + 1) % n
        p = gcd(A-B,n)
        if 1 < p < n: return p
        if p == n: return n

# Método de factorización de curva elíptica de Lenstra
def Lenstra(n: int, B: int = 1000, attempts: int = 100) -> int:
    for _ in range(attempts):
        
        x = randrange(n)
        y = randrange(n)
        a = randrange(n)
        
        try:
            
            b = (pow(y, 2, n) - (pow(x, 3, n) + (a * x) % n)) % n
            c = EllipticCurve(a, b, n)
            p = Point(c, x, y)
            
            
            current = p
            for i in range(2, B + 1):
                current = current * i  
                
        except ZeroDivisionError as e:
            
            den = int(str(e).split()[0])
            factor = gcd(den, n)
            if 1 < factor < n:
                return factor
                
    return None


def main():
    n = 455459
    factor = Lenstra(n)
    if factor:
        print(f"Un factor de n es: {factor}")
    else:
        print("No se encontró factor")

if __name__ == "__main__":
    main()
