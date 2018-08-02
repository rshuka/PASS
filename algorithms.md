# Klassifizierung der Optimierungsverfahren

Optimiserungsverfahren lassen sich zunächst einmal einteilen in exakte Optimierungsverfahren, die garantiert die Optimallösung finden, allerdings aus Komplexitätsgründen für schwierige Probleme nicht in Frage kommen, und Heuristiken, die nicht immer die Optimallösung finden. Die meisten Menschen denken bei dem Begriff Heuristik allerdings an die klassischen deterministischen Heuristiken, wie z.B. Prioritätsregeln von der Art „mache den dringendsten Job immer zuerst“. Naturanaloge Optimierungsverfahren sind zwar auch Heuristiken in dem Sinn, dass sie nicht immer die Optimallösung finden, aber sie funktionieren auf ganz anderen Prinzipien.

## Klassifizierung der Natuanalogen Optimierungsverfahren

- Trajektorenbasierte Optimierungsverfahren
  - Compass Search (CS)
  - Corana’s Simulated Annealing (SA)
- Populationsbasierte Optimierungsverfahren
  - Covariance Matrix Adaptation Strategy (CMA-ES)
  - Differential Evolution (DE)
  - Simple Genetic Algorithm
  - Artificial Bee Colony (ABC)

### Trajektorenbasierte Optimierungsverfahren

Trajektorenbasierte Optimierungsverfahren starten mit einem zufällig gewählten Agent und verfeinern die Lösung indem die aktuell gefundene Lösung mit einer besseren in der Nachbarschaft gefundene Lösung ersetzt wird.
