package com.gnormu.battleship.strategy;

import java.util.Arrays;
import com.gnormu.battleship.config.GameConfig;
import com.gnormu.battleship.domain.BoardView;
import com.gnormu.battleship.domain.CellContent;

public class MontecarloStrategy extends AbstractBattleshipStrategy {

    private static final byte DIM = GameConfig.BOARD_DIMENSION;

    private static final int[] INITIAL_HEAT_MAP = buildInitialHeatMap();

    /**
     * Construye el mapa de calor inicial evaluando las posiciones válidas para
     * todos los barcos sobre un tablero limpio.
     *
     * @return Mapa de calor con las posiciones válidas de los barcos
     */
    private static int[] buildInitialHeatMap() {
        int[] map = new int[GameConfig.DIMENSION_SQUARED];
        for (int cellIdx = 0; cellIdx < map.length; cellIdx++) {
            for (byte ship = 1; ship < CellContent.SHIP_SIZE; ship++) {
                map[cellIdx] += calculateInitialProbabilities(cellIdx, CellContent.getLength(ship));
            }
        }
        return map;
    }

    /**
     * Cuenta posiciones válidas (horizontal + vertical) para un barco de longitud
     * {@code shipLength} sobre la celda {@code cellIdx} en un tablero limpio.
     *
     * @param cellIdx    índice que representa la coordenada
     * @param shipLength largo del barco
     * @return Cantidad de posiciones válidas para el barco
     */
    private static int calculateInitialProbabilities(int cellIdx, int shipLength) {
        if (shipLength <= 0) {
            return 0;
        }
        int row = cellIdx / DIM;
        int col = cellIdx % DIM;
        int horizontal = Math.max(0, Math.min(col, DIM - shipLength) - Math.max(0, col - shipLength + 1) + 1);
        int vertical = Math.max(0, Math.min(row, DIM - shipLength) - Math.max(0, row - shipLength + 1) + 1);
        return horizontal + vertical;
    }

    private final int[] heatMap = new int[GameConfig.DIMENSION_SQUARED];
    private byte lastCoord = Byte.MIN_VALUE;

    private byte currentTargetShip = -1;
    private int currentHitCount = 0;
    private int minHitCell = -1;
    private int maxHitCell = -1;

    public MontecarloStrategy() {
        reset();
    }

    /** {@inheritDoc} */
    @Override
    public byte calculateNextShot(BoardView boardView) {
        if (lastShoot > 0) {
            onHit(lastCoord, (byte) lastShoot, boardView);
        } else if (lastCoord != Byte.MIN_VALUE) {
            onMiss(lastCoord, boardView);
        }

        applyHuntTargetBoost(boardView);

        lastCoord = getBestCoordFrom(heatMap);
        heatMap[lastCoord] = Byte.MIN_VALUE;
        return lastCoord;
    }

    /** {@inheritDoc} */
    @Override
    public void resetStrategy() {
        System.arraycopy(INITIAL_HEAT_MAP, 0, heatMap, 0, heatMap.length);
        lastCoord = Byte.MIN_VALUE;
        currentTargetShip = -1;
        currentHitCount = 0;
        minHitCell = -1;
        maxHitCell = -1;
    }

    private void onMiss(int coord, BoardView boardView) {
        int row = coord / DIM;
        int col = coord % DIM;

        if (currentTargetShip == -1) {
            for (byte ship = 1; ship < CellContent.SHIP_SIZE; ship++) {
                if (!boardView.isShipSunk(ship)) {
                    subtractInvalidatedWindows(coord, row, col, ship, 0, boardView);
                }
            }
        } else {
            subtractInvalidatedWindows(coord, row, col, currentTargetShip, currentHitCount, boardView);
        }
    }

    private void onHit(int coord, byte ship, BoardView boardView) {
        Arrays.fill(heatMap, 0);
        if (boardView.isShipSunk(ship)) {
            findNextTarget(boardView);
            if (currentTargetShip != -1) {
                scanAllWindowsForShip(currentTargetShip, boardView);
            } else {
                for (byte remaining = 1; remaining < CellContent.SHIP_SIZE; remaining++) {
                    if (!boardView.isShipSunk(remaining)) {
                        scanAllWindowsForShip(remaining, boardView);
                    }
                }
            }
        } else {
            if (currentTargetShip != ship) {
                recalculateTargetState(ship, boardView);
            } else {
                currentHitCount++;
                if (coord < minHitCell) {
                    minHitCell = coord;
                }
                if (coord > maxHitCell) {
                    maxHitCell = coord;
                }
            }
            scanAllWindowsForShip(currentTargetShip, boardView);
        }
    }

    private void subtractInvalidatedWindows(int coord, int row, int col,
            byte ship, int requiredHits, BoardView boardView) {
        int shipLength = CellContent.getLength(ship);
        byte hitState = (byte) -ship;

        for (int startCol = Math.max(0, col - shipLength + 1); startCol <= Math.min(col,
                DIM - shipLength); startCol++) {
            subtractWindowIfWasValid(row * DIM + startCol, 1, coord, shipLength, requiredHits, hitState, boardView);
        }
        for (int startRow = Math.max(0, row - shipLength + 1); startRow <= Math.min(row,
                DIM - shipLength); startRow++) {
            subtractWindowIfWasValid(startRow * DIM + col, DIM, coord, shipLength, requiredHits, hitState, boardView);
        }
    }

    private void subtractWindowIfWasValid(int windowStart, int stride, int changedCoord,
            int shipLength, int requiredHits, byte hitState, BoardView boardView) {
        int hitsFound = 0;
        for (int offset = 0; offset < shipLength; offset++) {
            int cellIdx = windowStart + offset * stride;
            if (cellIdx == changedCoord) {
                continue;
            }
            byte cellState = boardView.getCellState((byte) cellIdx);
            if (cellState == hitState) {
                hitsFound++;
            } else if (cellState != CellContent.WATER) {
                return;
            }
        }
        if (hitsFound != requiredHits) {
            return;
        }
        for (int offset = 0; offset < shipLength; offset++) {
            int cellIdx = windowStart + offset * stride;
            if (cellIdx != changedCoord && boardView.getCellState((byte) cellIdx) == CellContent.WATER) {
                heatMap[cellIdx]--;
            }
        }
    }

    private void scanAllWindowsForShip(byte ship, BoardView boardView) {
        int shipLength = CellContent.getLength(ship);
        byte hitState = (byte) -ship;

        if (ship != currentTargetShip || currentHitCount == 0) {
            scanHorizontalWindows(shipLength, hitState, 0, boardView);
            scanVerticalWindows(shipLength, hitState, 0, boardView);
        } else {
            int firstHitRow = minHitCell / DIM;
            int firstHitCol = minHitCell % DIM;
            for (int col = Math.max(0, firstHitCol - shipLength + 1); col <= Math.min(firstHitCol,
                    DIM - shipLength); col++) {
                scanWindow(firstHitRow * DIM + col, 1, currentHitCount, shipLength, hitState, boardView);
            }
            for (int row = Math.max(0, firstHitRow - shipLength + 1); row <= Math.min(firstHitRow,
                    DIM - shipLength); row++) {
                scanWindow(row * DIM + firstHitCol, DIM, currentHitCount, shipLength, hitState, boardView);
            }
        }
    }

    private void scanHorizontalWindows(int shipLength, byte hitState, int requiredHits, BoardView boardView) {
        for (int row = 0; row < DIM; row++) {
            int rowStart = row * DIM;
            for (int col = 0; col <= DIM - shipLength; col++) {
                scanWindow(rowStart + col, 1, requiredHits, shipLength, hitState, boardView);
            }
        }
    }

    private void scanVerticalWindows(int shipLength, byte hitState, int requiredHits, BoardView boardView) {
        for (int col = 0; col < DIM; col++) {
            for (int row = 0; row <= DIM - shipLength; row++) {
                scanWindow(row * DIM + col, DIM, requiredHits, shipLength, hitState, boardView);
            }
        }
    }

    private void recalculateTargetState(byte targetShip, BoardView boardView) {
        currentTargetShip = targetShip;
        currentHitCount = 0;
        minHitCell = -1;
        maxHitCell = -1;
        byte hitState = (byte) -targetShip;
        for (int cellIdx = 0; cellIdx < GameConfig.DIMENSION_SQUARED; cellIdx++) {
            if (boardView.getCellState((byte) cellIdx) == hitState) {
                if (currentHitCount == 0) {
                    minHitCell = cellIdx;
                }
                maxHitCell = cellIdx;
                currentHitCount++;
            }
        }
    }

    private void findNextTarget(BoardView boardView) {
        for (byte ship = 1; ship < CellContent.SHIP_SIZE; ship++) {
            if (!boardView.isShipSunk(ship)) {
                recalculateTargetState(ship, boardView);
                if (currentHitCount > 0) {
                    return;
                }
            }
        }
        currentTargetShip = -1;
        currentHitCount = 0;
        minHitCell = -1;
        maxHitCell = -1;
    }

    private void applyHuntTargetBoost(BoardView boardView) {
        if (currentTargetShip == -1) {
            return;
        }
        if (currentHitCount == 1) {
            updateNeighborsOfSingleHit(minHitCell, currentTargetShip, boardView);
        } else if (currentHitCount >= 2) {
            updateNeighborsOfMultipleHits(minHitCell, maxHitCell, currentTargetShip, boardView);
        }
    }

    private void updateNeighborsOfSingleHit(int hitCoord, byte targetShip, BoardView boardView) {
        int row = hitCoord / DIM;
        int col = hitCoord % DIM;
        int shipLength = CellContent.getLength(targetShip);
        byte hitState = (byte) -targetShip;

        int left = countWalkable(row * DIM + col - 1, -1, hitState, boardView);
        int right = countWalkable(row * DIM + col + 1, 1, hitState, boardView);
        if (left + right + 1 >= shipLength) {
            if (left > 0) {
                heatMap[hitCoord - 1] = Byte.MAX_VALUE;
            }
            if (right > 0) {
                heatMap[hitCoord + 1] = Byte.MAX_VALUE;
            }
        }

        int up = countWalkable((row - 1) * DIM + col, -DIM, hitState, boardView);
        int down = countWalkable((row + 1) * DIM + col, DIM, hitState, boardView);
        if (up + down + 1 >= shipLength) {
            if (up > 0) {
                heatMap[hitCoord - DIM] = Byte.MAX_VALUE;
            }
            if (down > 0) {
                heatMap[hitCoord + DIM] = Byte.MAX_VALUE;
            }
        }
    }

    private void updateNeighborsOfMultipleHits(int minHit, int maxHit, byte targetShip, BoardView boardView) {
        if ((minHit % DIM) == (maxHit % DIM)) {
            updateVerticalNeighbors(minHit, maxHit, targetShip, boardView);
        } else {
            updateHorizontalNeighbors(minHit, maxHit, targetShip, boardView);
        }
    }

    private void updateHorizontalNeighbors(int minHit, int maxHit, byte targetShip, BoardView boardView) {
        int row = minHit / DIM;
        int minCol = minHit % DIM;
        int maxCol = maxHit % DIM;
        byte hitState = (byte) -targetShip;
        int shipLength = CellContent.getLength(targetShip);

        int left = countWalkable(row * DIM + minCol - 1, -1, hitState, boardView);
        int right = countWalkable(row * DIM + maxCol + 1, 1, hitState, boardView);
        int hitSpan = maxCol - minCol + 1;
        if (left + right + hitSpan >= shipLength) {
            if (left > 0) {
                heatMap[row * DIM + (minCol - 1)] = Byte.MAX_VALUE;
            }
            if (right > 0) {
                heatMap[row * DIM + (maxCol + 1)] = Byte.MAX_VALUE;
            }
        }
    }

    private void updateVerticalNeighbors(int minHit, int maxHit, byte targetShip, BoardView boardView) {
        int col = minHit % DIM;
        int minRow = minHit / DIM;
        int maxRow = maxHit / DIM;
        byte hitState = (byte) -targetShip;
        int shipLength = CellContent.getLength(targetShip);

        int up = countWalkable((minRow - 1) * DIM + col, -DIM, hitState, boardView);
        int down = countWalkable((maxRow + 1) * DIM + col, DIM, hitState, boardView);
        int hitSpan = maxRow - minRow + 1;
        if (up + down + hitSpan >= shipLength) {
            if (up > 0) {
                heatMap[(minRow - 1) * DIM + col] = Byte.MAX_VALUE;
            }
            if (down > 0) {
                heatMap[(maxRow + 1) * DIM + col] = Byte.MAX_VALUE;
            }
        }
    }

    /**
     * Cuenta celdas transitables (agua o impacto del barco objetivo) desde
     * {@code startIdx} en dirección {@code stride} hasta salir del tablero o
     * encontrar una celda bloqueada.
     */
    private int countWalkable(int startIdx, int stride, byte hitState, BoardView boardView) {
        int count = 0;
        for (int cellIdx = startIdx; cellIdx >= 0 && cellIdx < GameConfig.DIMENSION_SQUARED; cellIdx += stride) {
            byte cellState = boardView.getCellState((byte) cellIdx);
            if (cellState == CellContent.WATER || cellState == hitState) {
                count++;
            } else {
                break;
            }
        }
        return count;
    }

    /**
     * Evalúa una ventana de {@code shipLength} celdas desde {@code windowStart} con
     * paso {@code stride}. Si contiene exactamente {@code requiredHits} impactos y
     * ningún bloqueo, incrementa el heatmap para las celdas de agua de la ventana.
     */
    private void scanWindow(int windowStart, int stride, int requiredHits, int shipLength,
            byte hitState, BoardView boardView) {
        boolean valid = true;
        int hitsFound = 0;
        for (int offset = 0; offset < shipLength; offset++) {
            byte cellState = boardView.getCellState((byte) (windowStart + offset * stride));
            if (cellState == hitState) {
                hitsFound++;
            } else if (cellState != CellContent.WATER) {
                valid = false;
                break;
            }
        }
        if (valid && hitsFound == requiredHits) {
            for (int offset = 0; offset < shipLength; offset++) {
                int cellIdx = windowStart + offset * stride;
                if (boardView.getCellState((byte) cellIdx) == CellContent.WATER) {
                    heatMap[cellIdx]++;
                }
            }
        }
    }

    private byte getBestCoordFrom(int[] currentHeatMap) {
        int maxProbability = Integer.MIN_VALUE;
        byte bestCoord = Byte.MIN_VALUE;
        for (int cellIdx = 0; cellIdx < GameConfig.DIMENSION_SQUARED; cellIdx++) {
            if (currentHeatMap[cellIdx] > maxProbability) {
                maxProbability = currentHeatMap[cellIdx];
                bestCoord = (byte) cellIdx;
            }
        }
        return bestCoord;
    }
}