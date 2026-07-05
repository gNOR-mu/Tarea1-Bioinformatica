package com.gnormu.battleship.strategy;

import java.util.Arrays;
import com.gnormu.battleship.config.GameConfig;
import com.gnormu.battleship.domain.BoardView;
import com.gnormu.battleship.domain.CellContent;

public class Montecarlo extends AbstractBattleshipStrategy {

    private static final byte DIM = GameConfig.BOARD_DIMENSION;

    private static final int[] INITIAL_HEAT_MAP = buildInitialHeatMap();

    private static int[] buildInitialHeatMap() {
        int[] map = new int[GameConfig.DIMENSION_SQUARED];
        for (int i = 0; i < map.length; i++) {
            for (byte ship = 1; ship < CellContent.SHIP_SIZE; ship++) {
                map[i] += calculateInitialProbabilities(i, CellContent.getLength(ship));
            }
        }
        return map;
    }

    private static int calculateInitialProbabilities(int idx, int length) {
        if (length <= 0)
            return 0;
        int row = idx / DIM;
        int col = idx % DIM;

        int minCol = Math.max(0, col - length + 1);
        int maxCol = Math.min(col, DIM - length);
        int horizontal = Math.max(0, maxCol - minCol + 1);

        int minRow = Math.max(0, row - length + 1);
        int maxRow = Math.min(row, DIM - length);
        int vertical = Math.max(0, maxRow - minRow + 1);

        return horizontal + vertical;
    }

    private final int[] heatMap = new int[GameConfig.DIMENSION_SQUARED];
    private final byte[] boardMemory = new byte[GameConfig.DIMENSION_SQUARED];
    private final byte[] hitCounts = new byte[CellContent.SHIP_SIZE];
    private boolean isFirstTurn;
    private byte lastCoord;

    public Montecarlo() {
        reset();
    }

    @Override
    public byte calculateNextShot(BoardView boardView) {
        if (lastCoord != -1) {
            if (lastShoot > 0) {
                boardMemory[lastCoord] = (byte) -lastShoot;
                hitCounts[lastShoot]++;
            } else {
                boardMemory[lastCoord] = CellContent.MISS;
            }
        }

        byte nextShot;
        if (isFirstTurn) {
            nextShot = getBestCoordFrom(INITIAL_HEAT_MAP);
        } else {
            recalculateHeatMap(boardView);
            nextShot = getBestCoordFrom(heatMap);
        }

        lastCoord = nextShot;
        return nextShot;
    }

    @Override
    public void resetStrategy() {
        System.arraycopy(INITIAL_HEAT_MAP, 0, heatMap, 0, heatMap.length);
        Arrays.fill(boardMemory, CellContent.WATER);
        Arrays.fill(hitCounts, (byte) 0);
        isFirstTurn = true;
        lastCoord = -1;
    }

    private void recalculateHeatMap(BoardView boardView) {
        Arrays.fill(heatMap, 0);

        byte targetShip = -1;
        for (byte ship = 1; ship < CellContent.SHIP_SIZE; ship++) {
            if (hitCounts[ship] > 0 && !boardView.isShipSunk(ship)) {
                targetShip = ship;
                break;
            }
        }

        if (targetShip != -1) {
            calculateProbabilitiesForShip(targetShip, hitCounts[targetShip]);
        } else {
            for (byte ship = 1; ship < CellContent.SHIP_SIZE; ship++) {
                if (!boardView.isShipSunk(ship)) {
                    calculateProbabilitiesForShip(ship, 0);
                }
            }
        }
    }

    private void calculateProbabilitiesForShip(byte ship, int requiredHits) {
        int length = CellContent.getLength(ship);
        byte targetState = (byte) -ship;

        for (int row = 0; row < DIM; row++) {
            int rowOffset = row * DIM;

            for (int col = 0; col <= DIM - length; col++) {
                boolean valid = true;
                int hitsFound = 0;

                for (int k = 0; k < length; k++) {
                    byte state = boardMemory[rowOffset + col + k];

                    if (state == targetState) {
                        hitsFound++;
                    } else if (state != CellContent.WATER) {
                        valid = false;
                        break;
                    }
                }

                if (valid && hitsFound == requiredHits) {
                    for (int k = 0; k < length; k++) {
                        int index = rowOffset + col + k;
                        if (boardMemory[index] == CellContent.WATER) {
                            heatMap[index]++;
                        }
                    }
                }
            }
        }

        for (int col = 0; col < DIM; col++) {
            for (int row = 0; row <= DIM - length; row++) {
                boolean valid = true;
                int hitsFound = 0;

                for (int k = 0; k < length; k++) {
                    byte state = boardMemory[(row + k) * DIM + col];

                    if (state == targetState) {
                        hitsFound++;
                    } else if (state != CellContent.WATER) {
                        valid = false;
                        break;
                    }
                }

                if (valid && hitsFound == requiredHits) {
                    for (int k = 0; k < length; k++) {
                        int index = (row + k) * DIM + col;
                        if (boardMemory[index] == CellContent.WATER) {
                            heatMap[index]++;
                        }
                    }
                }
            }
        }
    }

    /**
     * Extrae la coordenada con la mayor densidad de probabilidad.
     */
    private byte getBestCoordFrom(int[] currentHeatMap) {
        isFirstTurn = false;
        int maxProbability = -1;
        byte bestCoord = -1;

        for (int i = 0; i < GameConfig.DIMENSION_SQUARED; i++) {
            if (currentHeatMap[i] > maxProbability) {
                maxProbability = currentHeatMap[i];
                bestCoord = (byte) i;
            }
        }

        return bestCoord;
    }
}