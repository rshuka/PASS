# Artificial Bee Colony

pagmo2 implementiert exakt die Version des artifical bee colony Optimierers aus [https://abc.erciyes.edu.tr/pub/NevImpOfABC.pdf](), Algorithm 2. ABC ist ein populationsbasiertes Optimierungsverfahren für single-objective unconstrained Optimierungsprobleme.

Der Algorithmus kann über zwei Parameter konfiguriert werden: Die Populationsgröße (`SN`), und die Anzahl an Fehlversuchen (`limit`), bevor ein Agent von "beschäftigt" auf "erkunden" wechselt. Er arbeitet nach folgendem Schema:

    Erzeuge eine Population mit SN Agenten
    Für jeden Agenten:
        Initialisiere Agent
    Bis maximale Anzahl Evaluationen erreicht ist:
        Falls Agenten existieren mit Fehlversuchszähler ≥ limit:
            Wähle einen einzigen Agenten mit größtem Fehlversuchszähler als Kundschafter
        Für jeden Agenten außer dem Kundschafter:
            // 1. Employed Bees Phase
            Aktualisiere Agent
        Speichere die beste bisher gefundene Lösung
        Falls ein Kundschafter existiert:
            // 2. Scout Bees Phase
            Re-initialisiere den Kundschafter
        Wiederhole SN mal:
            // 3. Onlooker Bees Phase
            Wähle einen zufälligen Agenten m mit Wahrscheinlichkeit p_m
            Aktualisiere Agent m
        Abspeichern der besten bisher gefundenen Lösung


## Initialisierung

Wird ein Agent (re-)initialisiert, wird ihm eine gleichverteilt zufällige Position innerhalb der Problemgrenzen zugewiesen, und sein Fehlversuchzähler wird auf 0 gesetzt.

## Aktualisierung eines Agenten

Die Aktualisierung eines Agenten läuft in den folgenden Schritten ab:
1. Sei die aktuelle Position des Agenten `x_m`.
2. Wähle eine zufällige Dimension _d_.
3. Wähle einen zufälligen anderen Agenten; seine Position sei `x_k`.
4. Wähle ein zufälliges ϕ∈[-1, 1].
5. Berechne eine neue Position `v_m`, indem ϕ * (`x_m_d` - `x_k_d`) zum Wert in Dimension _d_ der aktuellen Position `x_m` addiert wird.
6. Ist f(`v_m`) < f(`x_m`), bewege den Agenten nach `v_m` und setze seinen Fehlversuchszähler auf 0; sonst erhöhe seinen Fehlversuchszähler um 1.

## Wahl eines zufälligen Agenten

Die Wahrscheinlichkeit `p_m`, mit der Agent `x_m` während der _onlooker bees_ Phase ausgewählt wird, beträgt

                     ⎛ SN         ⎞-1
    p_m = fit(x_m) * ⎜ ∑  fit(x_i)⎟
                     ⎝i=1         ⎠

wobei

               ⎧ 1 / (1+f(x_m)) , falls f(x_m) ≥ 0
    fit(x_m) = ⎨
               ⎩ 1 + |f(x_m)|   , falls f(x_m) < 0

Dieser Fitnesswert hat im Gegensatz zum objective value `f` die nützliche Eigenschaft, dass stets `fit(x_m)` ≥ 1 gilt.

## Zusammenfassung

Die erste Phase _employed bees_ weist Ähnlichkeiten zum PSO auf: Jeder Agent bewegt sich auf genau einen zufälligen anderen Agenten zu. Anders als beim PSO bewegen sich ABC Agenten jedoch nur in einer einzigen Dimension; außerdem bleiben sie stehen, anstatt zu einer Position mit einem schlechteren objective value zu wechseln.

Die zweite Phase _scout bees_ implementiert ein Neustart-Kriterium, das dann einsetzt, wenn sich ein Agent für mindestens `limit` Evaluationen nicht bewegt hat. Damit soll der _exploration_-Aspekt der Optimierung gewährleistet werden. Pro Iteration kann nur ein einzier Agent neu gestartet werden.

Die dritte Phase _onlooker bees_ dient dazu, eine gute _exploitation_ zu erreichen, indem diejenigen Agenten die meisten Evaluationen erhalten, die global betrachtet an den besten bisher gefundenen Positionen stehen. Da in dieser Phase noch einmal `SN` Evaluationen durchgeführt werden, ist die Population quasi 2 * `SN` groß, wobei nur die Hälfte dieser Population mit einem Gedächtnis versehen ist.
