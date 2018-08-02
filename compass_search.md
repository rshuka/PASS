# Compass Search

Der Algorithmus wurde 2003 vorgestellt [https://epubs.siam.org/doi/pdf/10.1137/S003614450242889]. Compass Search ist einer der ersten Beispielen weche auf Direct Search [Cit Hooke Jeeves] basieren und auf einem Computer implementiert wurden um einem Optimierungsproblem zu lösen.

## Initialisierung

1. `k` ist der Index für die Iterationen
2. `x_k ∈ R^n` is das k-te element
3. `x_0` ist der erster Agent
4. `D⊕` ist das Koordinaten System. Die Vektoren sind positiv und negativ.
    `D⊕ = {e_1,e_2,...,e_n,−e_1,−e_2,...,−e_n}`
5. `∆_k` Schrittlänge Parameter welcher die Länge des Schrittes kontrolliert. Bei der Initialisierung muss `∆_0` definiert werden. Als empfohlener Wert wird `∆_0 = 0.3` genommen. Bedingung: `∆_0` muss größer als der Toleranzwert sein.
6. `∆_tol > 0` Ein Toleranzwert für die Konvergenz des Problems

## Pseudo Code

Für jede Iteration `K = 1, 2, ...` wiederhole die Schritte bis das Optimum gefunden wird.
```{r, tidy=FALSE, eval=FALSE, highlight=FALSE }

Schritt 1: Lass `D⊕` den Satz von Koordinatenrichtungen sein {±e_i |i = 1,...,n}, wo e_i die i-te Koordinate in R^n ist

Schritt 2: Wenn eine d_k ∈ D⊕ existiert welche f(x_k + ∆_kd_k) < f(x_k), dann mach folgendes    -> wenn ein neues besseres Lösung gefunden wird
          - Setze x_k+1 = x_k + ∆_kd_k (Iteration wechseln)
          - Setze ∆_k+1 = ∆_k (keine Änderung des Schrittlängen Parameters)

Schritt 3: Sonst,f(x_k + ∆_kd) ≥ f(x_k) für alle d ∈ D⊕, dann nache folgendes                   -> wenn keine bessere Lösung gefunden wird
          - Setze x_k+1 = x_k (keine Iteration wechseln)
          - Setze ∆_k+1 = 1/2 * ∆_k (Änderung des Schrittlängen Parameters)
          - Wenn ∆_k+1 < ∆_tol, dann abbrechen.

```

## Zusammenfassung

Ganz am Anfang wird ein zufälliges Punk in den Suchraum gewählt. Als erstes Schritt werden alle Nachbarn in jeder Richtung bestimmt. Dieser Nachbarn werden geprüft und wenn einer der Nachbarn eine bessere Lösung liefert, wird dieser gewählt und von vorne angefangen. Wenn eine bessere Lösung gefunden wird, wird die Schrittlänge nicht geändert.
Wird keine bessere Lösung gefunden, dann wird die Schrittlänge halbiert und nochmal alle Nachbarn geprüfuft.
Das passiert bis die Schrittlänge einem vordefinierten Minimum erreicht hat. Dieser Minimum muss kleiner als der Anfangswert sein.
Wenn die gewünsche Lösung schon vorher gefunden wird, wird auch abgebrochen.
