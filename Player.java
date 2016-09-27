package slather.g4;

import slather.sim.Cell;
import slather.sim.Point;
import slather.sim.Move;
import slather.sim.Pherome;
import java.util.*;

public class Player implements slather.sim.Player {

	private Random gen;
	private final static int ANGEL_RANGE = 360;
	private final static int NUMBER_OF_RANDOM_TRY = 4;
	private final static double PHEROME_IMPORTANCE = 0.2;
	private int tail;
	private double visible_distance;

	public void init(double d, int t, int side_length) {
		this.gen = new Random();
		this.visible_distance = d;
		this.tail = t;
	}

	public Move play(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		if (player_cell.getDiameter() >= 2) // reproduce whenever possible
			return new Move(true, (byte) -1, (byte) -1);

		// care bout your children, dude!

		// strategies choosen branch
		int arg = 0;
		int spin_sep = get_spin_sep(this.tail, this.visible_distance); 
		int detector_sep = get_detector_sep(this.tail, this.visible_distance);
		if (nearby_cells.size() == 0) {
			arg = memory;
		} else if (isCrowded(player_cell, nearby_cells, nearby_pheromes) == true) {
			arg = spin(player_cell, memory, nearby_cells, nearby_pheromes, spin_sep);
		} else {
			arg = detector(player_cell, memory, nearby_cells, nearby_pheromes, detector_sep);
		}

		Point vector = extractVectorFromAngle(arg);
		if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) {
			return new Move(vector, (byte) arg);
		}

		// backup strategy: random escape
		for (int i = 0; i < Player.NUMBER_OF_RANDOM_TRY; i++) {
			arg = gen.nextInt(Player.ANGEL_RANGE) + 1;
			vector = extractVectorFromAngle(arg);
			if (!collides(player_cell, vector, nearby_cells, nearby_pheromes))
				return new Move(vector, (byte) arg);
		}
		return new Move(new Point(0, 0), (byte) 0);	
	}

	private int get_spin_sep(int tail, double visible_distance) {
		return 10;
	}

	private int get_detector_sep(int tail, double visible_distance) {
		return 4;
	}

	private double threshold(int tail, double visible_distance) {
		return 2;
	}

	private int spin(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes, int seperation) {
		return memory + Player.ANGEL_RANGE / seperation;
	}

	private int detector(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes, int seperation) {
		double[] direction = new double[seperation];

		for (Cell nearby_cell : nearby_cells) {
			double m_x = nearby_cell.getPosition().x - player_cell.getPosition().x;
			double m_y = nearby_cell.getPosition().y - player_cell.getPosition().y;
			double theta = Math.atan(m_y/m_x);
			int index = (int) ((theta / (2 * Math.PI) + (m_x > 0 ? (m_y < 0 ? 1 : 0) : 0.5)) * seperation);
			direction[index] += trans1(nearby_cell.distance(player_cell));
		}

		for (Pherome nearby_pherome : nearby_pheromes) {
			if (nearby_pherome.player == player_cell.player) continue;
			double m_x = nearby_pherome.getPosition().x - player_cell.getPosition().x;
			double m_y = nearby_pherome.getPosition().y - player_cell.getPosition().y;
			double theta = Math.atan(m_y/m_x);
			int index = (int) ((theta / (2 * Math.PI) + (m_x > 0 ? (m_y < 0 ? 1 : 0) : 0.5)) * seperation);
			direction[index] += Player.PHEROME_IMPORTANCE * trans1(nearby_pherome.distance(player_cell));
		}

		double sum = 0;

		for (int index = 0; index < direction.length; ++index) {
			sum += trans2(direction[index]);
		}

		// We give each direction a probability direction[i]/sum
		sum *= Math.random();
		int index = 0;
		for (; index < direction.length; ++index) {
			sum -= trans2(direction[index]);
			if (sum <= 0) break;
		}

		int arg = 0;
		return arg = index * Player.ANGEL_RANGE / seperation + gen.nextInt(Player.ANGEL_RANGE / seperation);
	}

	private boolean isCrowded(Cell player_cell, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		double weightSum = 0;
		for (Cell nearby_cell : nearby_cells) {
			weightSum += trans1(nearby_cell.distance(player_cell));
		}

		for (Pherome nearby_pherome : nearby_pheromes) {
			if (nearby_pherome.player != player_cell.player) {
				weightSum += Player.PHEROME_IMPORTANCE * trans1(nearby_pherome.distance(player_cell));
			}
		}

		return weightSum >= threshold(this.tail, this.visible_distance);
	}

	private double trans1(double a) {
		return 1.0/(a+0.5);
	}

	private double trans2(double a) {
		return 1/Math.pow((0.01+a),2);
	}

	private boolean collides(Cell player_cell, Point vector, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		Iterator<Cell> cell_it = nearby_cells.iterator();
		Point destination = player_cell.getPosition().move(vector);
		while (cell_it.hasNext()) {
			Cell other = cell_it.next();
			if (destination.distance(other.getPosition()) < 0.5 * player_cell.getDiameter() + 0.5 * other.getDiameter()
					+ 0.00011)
				return true;
		}
		Iterator<Pherome> pherome_it = nearby_pheromes.iterator();
		while (pherome_it.hasNext()) {
			Pherome other = pherome_it.next();
			if (other.player != player_cell.player
					&& destination.distance(other.getPosition()) < 0.5 * player_cell.getDiameter() + 0.0001)
				return true;
		}
		return false;
	}

	private Point extractVectorFromAngle(int arg) {
		double theta = Math.toRadians((double) arg);
		double dx = Cell.move_dist * Math.cos(theta);
		double dy = Cell.move_dist * Math.sin(theta);
		return new Point(dx, dy);
	}

}
