import numpy as np
import matplotlib.pyplot as plt


def oblicz_dEdx(energia, Z1, Z2, I, n_m3, e, me, c, mp):
    """
    Funkcja obliczająca wartość dE/dx (w J/m) dla protonu poruszającego się w materiale,
    zgodnie z uproszczoną wersją równania Bethe-Blocha.

    Parametry
    ---------
    energia : float
        Energia cząstki w danym momencie (w dżulach, SI).
    Z1 : int
        Liczba atomowa (ładunek) cząstki przychodzącej (dla protonu Z1=1).
    Z2 : int
        Liczba atomowa materiału docelowego (dla miedzi Z2=29).
    I : float
        Średni potencjał jonizacji (w dżulach, SI).
    n_m3 : float
        Gęstość elektronów w elektronach na metr sześcienny (elektrony/m³).
    e : float
        Ładunek elementarny (w kulombach, C).
    me : float
        Masa elektronu (w kg).
    c : float
        Prędkość światła (w m/s).
    mp : float
        Masa protonu (w kg).

    Zwraca
    -------
    dEdx : float
        Strata energii na jednostkę drogi (J/m). Formalnie powinna być ujemna,
        ale tutaj zwracamy wartość dodatnią, aby łatwo odjąć ją od energii w pętli.
    """
    # Zabezpieczenie na wypadek energii ujemnej lub zerowej (brak sensu fizycznego)
    if energia <= 0:
        return 0.0

    # Obliczenie prędkości (uproszczenie: formuła nierelatywistyczna v = sqrt(2E/m))
    v = np.sqrt(2 * energia / mp)  # m/s
    beta = v / c  # beta = v/c

    # Składowa logarytmiczna (argument musi być bezwymiarowy)
    ln_term = np.log((2 * me * v ** 2) / I)

    # Uproszczone równanie Bethe-Blocha:
    # dE/dx = (4π * Z1^2 * e^4 / (me * v^2)) * (n / Z2) * [ln(...) - ln(1 - beta^2) - beta^2]
    dEdx = (4 * np.pi * Z1 ** 2 * e ** 4 / (me * v ** 2)) \
           * (n_m3 / Z2) \
           * (ln_term - np.log(1 - beta ** 2) - beta ** 2)

    # Jeżeli dEdx wyszłoby ujemne (np. w rejonie niskich energii), zwróć zero
    if dEdx < 0:
        return 0.0

    return dEdx


# ------------------------------
# Główna część symulacji

# Parametry cząstki i materiału
Z1 = 1  # Proton
Z2 = 29  # Miedź
I = 322 * 1.602e-19  # Konwersja 322 eV -> dżule

# Gęstość elektronów: 8,5e22 elektronów na cm^3 => przelicznik: 1 cm^3 = 1e-6 m^3
n_cm3 = 8.5e22
n_m3 = n_cm3 * 1e6  # elektronów / m^3

e = 1.602e-19  # Ładunek elementarny w kulombach (C)
me = 9.109e-31  # Masa elektronu (kg)
mp = 1.67e-27  # Masa protonu (kg)
c = 3.0e8  # Prędkość światła (m/s)

# Energia początkowa ~3 MeV w dżulach (1 eV ~ 1.602e-19 J => 3 MeV ~ 4.8e-13 J)
energia_startowa = 5e-13
energia = energia_startowa

# Krok obliczeniowy w metrach (0.1 cm = 1.0e-3 m)
krok_x_m = 1e-3

# Minimalna energia (prog zatrzymania)
energia_min = 0.0
# Początkowa pozycja w metrach
x_m = 0.0

# Listy do zapisu wyników
pozycje_m = [x_m]
energie = [energia]

index = 0
limit = 1000000

# Główna pętla obliczeniowa
while energia > energia_min:
    index += 1
    # Obliczenie dE/dx (J/m)
    dEdx = oblicz_dEdx(energia, Z1, Z2, I, n_m3, e, me, c, mp)

    # Zabezpieczenie, jeśli dEdx jest zerowe (brak dalszej straty energii)
    if dEdx == 0.0:
        break

    # Strata energii na odcinku 'krok_x_m'
    deltaE = dEdx * krok_x_m
    energia_nowa = energia - deltaE

    # Unikamy zejścia poniżej zera (energia nie może być ujemna)
    if energia_nowa < 0.0:
        energia_nowa = 0.0

    # Uaktualnienie energii i pozycji
    energia = energia_nowa
    x_m += krok_x_m

    # Zapis wyników
    pozycje_m.append(x_m)
    energie.append(energia)
    if index == limit:
        break

# Ostateczna głębokość zatrzymania w metrach, konwersja na cm
glebia_zatrzymania_cm = x_m * 100

print(f"Głębokość zatrzymania cząstki: {glebia_zatrzymania_cm:.2f} cm")
# Wizualizacja
plt.figure(figsize=(10, 6))
pozycje_cm = [pos * 100 for pos in pozycje_m]  # Konwersja osi X na cm
plt.plot(pozycje_cm, energie, label="Energia cząstki (J)")
plt.xlabel("Głębokość w materiale (cm)")
plt.ylabel("Energia (J)")
plt.title("Utrata energii cząstki naładowanej w materiale (Bethe-Bloch)")
plt.grid()
plt.legend()
plt.show()