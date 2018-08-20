# Simulated Annealing

Die in pagmo2 implementierte Version [1] (von 2017) von Simulated Annealing basiert auf einer Veröffentlichung von Corana [2] von 1987.

## Initialisierung

Diese Implementierung von SA in [1] kann über folgende Parameter konfiguriert werden:

Name|Beschreibung|Vorgeschlagener Wert in [1]|Vorgeschlagener Wert in [2]
----|------------|---------------------------|---------------------------
_Ts_|Starttemperatur|10|-
_Tf_|finale Temperatur|0.1|-
_N\_A_|Die Temperatur _T_ wird in dieser Anzahl Schritten von _Ts_ auf _Tf_ abgesenkt.|10|-
_N\_T_|Die Anzahl an Suchintervalländerungen pro Temperaturschritt _N\_A_.|1|max(100, 5 * n)
_N\_S_|Anzahl an Evaluationen pro _N\_T_.|20|20
_Vs_|Initiale (sowie maximale) Distanz, die sich die Suchposition in einem Schritt bewegen darf.|1|-

## Pseudocode

```
Sei f : ℝ^n -> ℝ die zu minimierende Funktion, lb ∈ ℝ^n die unteren Grenzen, rb ∈ ℝ^n die oberen Grenzen des Definitionsbereichs.
Sei x_opt ∈ ℝ^n ein zufälliger Punkt innerhalb der Funktionsgrenzen
Sei v ∈ ℝ^n = (Vs, ..., Vs)
// [1] erwartet als Parameter N_A und berechnet daraus r_T; [2] erwartet als Parameter r_T und benutzt N_A überhaupt nicht, sondern prüft ein benutzerdefiniertes Abbruchkriterium
Setze r_T = (Tf / Ts)^(1.0 / N_A)
Wiederhole N_A mal:
  x = x_opt
  Wiederhole N_T mal:
    n_d ∈ ℝ^n = (0, ..., 0)
    Wiederhole N_S mal:
      Wiederhole für jede Dimension d = 1, ..., n:
        // [2] benutzt hier eine Schleife und wirft so lange neue Zufallszahlen, bis diese in den Grenzen lb_d und ub_d liegen.
        Wähle eine gleichverteilte Zufallszahl r aus Intervall [max(x_d - v_d, lb_d), min(x_d + v_d, ub_d)]
        x' = (x_1, ..., x_{h-1}, r, x_{h+1}, ..., x_n)
        Wenn f(x') <= f(x_opt):
          x_opt = x'
        Wenn f(x') <= f(x), oder mit Wahrscheinlichkeit p = exp(-|f(x) - f(x')| / T):
          x = x'
          n_d = n_d + 1
    Wiederhole für jede Dimension d = 1, ..., n:
      // [1] setzt c_d konstant auf 2 für alle d; dieser Wert wird auch von [2] empfohlen.
            ⎧v_d * (1 + c_d * ((n_d/N_S - 0.6) / 0.4)) , wenn n_d > 0.6N_S
      v_d = ⎨v_d / (1 + c_d * ((0.4 - n_d/N_S) / 0.4)) , wenn n_d < 0.4N_S
            ⎩v_d                                       , sonst
      Wenn v_d > Vs:
        // In [2] nicht vorgesehen.
        v_d = Vs
  T = r_T * T
```

## Zusammenfassung

SA ist ein nicht populationsbasierter Optimierungsalgorithmus, der von einem einzigen Punkt ausgehend diesen Punkt jeweils in einer einzigen Dimension verschiebt. Hat dieser Punkt einen besseren Funktionswert, wird dieser als neuer Ausgangspunkt übernommen. Hat der neue Punkt einen schlechteren Funktionswert, wird er dennoch mit einer gewissen Wahrscheinlichkeit angenommen.

Zentrales Element des Algorithmus ist die aktuelle Temperatur _T_. Je größer der Wert dieser Variablen ist, um so höher ist die Wahrscheinlichkeit, dass ein solcher uphill move durchgeführt wird. Die Temperatur startet auf einem hohen Wert und wird schrittweise abgesenkt. Das soll verhindern, dass der Optimierer in einem lokalen Minimum stecken bleibt.

Ein weiterer wichtiger Aspekt ist die dimensionsweise separat abgespeicherte Schrittweite _v_. Neue Punkte werden in einer maximalen Reichweite von _v_ zum aktuellen Punkt gesucht. _v_ wird regelmäßig so angepasst, dass etwa 50% der Funktionsevaluationen akzeptiert werden (wobei sowohl bessere Funktionswerte als auch uphill moves mitgezählt werden).

Insgesamt arbeitet SA folgendermaßen: Der Algorithmus sucht ausgehend von der aktuellen Position neue Punkte. Alle _N\_S_ Evaluationen wird _v_ neu berechnet. Nachdem _v_ _N\_T_ mal aktualisiert wurde, wird die Temperatur abgesenkt. Nach _N\_A_ Temperaturabsenkungen terminiert der Algorithmus.

[1]: https://esa.github.io/pagmo2/docs/cpp/algorithms/simulated_annealing.html
[2]: http://people.sc.fsu.edu/~inavon/5420a/corana.pdf