package slather.g4;

import slather.sim.Cell;
import slather.sim.Point;
import slather.sim.Move;
import slather.sim.Pherome;
import java.util.*;


public class Player implements slather.sim.Player {
    
    private Random gen;

    public void init(double d, int t) {
	gen = new Random();
    }

    public Move play(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		if (player_cell.getDiameter() >= 2) // reproduce whenever possible
	    	return new Move(true, (byte)-1, (byte)-2);
	    /*
	    if (memory == -1) {
	    	return new Move(extractVectorFromAngle(0), (byte)-3);
	    } else if (memory == -2) {
	    	return new Move(extractVectorFromAngle(180), (byte)-4);
	    }

	    if (memory == -3) {
	    	return new Move(extractVectorFromAngle(0), (byte)0);
	    } else if (memory == -4) {
	    	return new Move(extractVectorFromAngle(180), (byte)180);
	    }*/
/*
	    if (nearby_cells.size() == 1) {
	    		for(Cell c : nearby_cells) {
	    			if (c.player == player_cell.player) {
	    				//System.out.println("1---------------------------------------");
	    				double newX = player_cell.getPosition().x - c.getPosition().x;
	    				double newY = player_cell.getPosition().y - c.getPosition().y;
	    				Point pos = new Point(newX,  newY);
	    				int arg = (int) (Math.atan2(newY, newX)/3.1415926*180) + (newX < 0 ? 180 : 0);
	    				//System.out.println(arg + "-----------------------------------------");
	    				Point vector = extractVectorFromAngle(arg);
	    				
	    				if (!collides( player_cell, vector, nearby_cells, nearby_pheromes)) {
	    					
	    					return new Move(vector, (byte) arg);
	    				}
						
	    			}
	    		}
	    		
	    } 
*/    
//System.out.println("3");
	    if (nearby_cells.size() == 0) {
	    	Point vector = extractVectorFromAngle( (int)memory);
	    		// check for collisions
	    	if (!collides( player_cell, vector, nearby_cells, nearby_pheromes))
				return new Move(vector, memory);
	    }

	    int[] direction = new int[4];
	    int min_dir = Integer.MAX_VALUE;
	    List<Integer> min_index = new ArrayList<>();

	    for (Cell nearby_cell : nearby_cells) {
	    	double m_x = nearby_cell.getPosition().x - player_cell.getPosition().x;
	    	double m_y = nearby_cell.getPosition().y - player_cell.getPosition().y;
	    	if (m_x >= 0 && m_y >= 0) { ++direction[0]; }
	    	if (m_x >= 0 && m_y < 0) { ++direction[3]; }
	    	if (m_x < 0 && m_y >= 0) { ++direction[1]; }
	    	if (m_x < 0 && m_y < 0) { ++direction[2]; }
	    	//System.out.println("distance: " + nearby_cell.distance(player_cell));
	    	//System.out.println("Diameter: " + nearby_cell.getDiameter());
	    }

	    for (int i = 0 ; i < direction.length ; ++i) {
	    	if (direction[i] <= min_dir) {
	    		min_dir = direction[i];
	    	}
	    }

	    for (int i = 0 ; i < direction.length ; ++i) {
	    	if (direction[i] == min_dir) {
	    		min_index.add(i);
	    	}
	    }

	    int selected = gen.nextInt(min_index.size());

	    //int arg = gen.nextInt(30) + min_index * 90 + 30;
	    int arg = min_index.get(selected) * 90 + gen.nextInt(90);

	    Point vector = extractVectorFromAngle(arg);
	   	if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) {
			return new Move(vector, (byte) arg);
		}
		for (int i=0; i<4; i++) {
	    	int arg2 = gen.nextInt(360)+1;
	    	Point vector2 = extractVectorFromAngle(arg2);
	    	if (!collides(player_cell, vector2, nearby_cells, nearby_pheromes)) 
				return new Move(vector2, (byte) arg2);
		}
		//System.out.println("4");
		return new Move(new Point(0,0), (byte)0);
	    /*
	    for (Cell nearby_p : nearby_pheromes) {
	    	double m_x = nearby_p.getPosition().x - player_cell.getPosition().x;
	    	double m_y = nearby_p.getPosition().y - player_cell.getPosition().y;
	    	if (m_x >= 0 && m_y >= 0) { ++direction[0]; }
	    	if (m_x >= 0 && m_y < 0) { ++direction[3]; }
	    	if (m_x < 0 && m_y >= 0) { ++direction[1]; }
	    	if (m_x < 0 && m_y < 0) { ++direction[2]; }
	    	//System.out.println("distance: " + nearby_cell.distance(player_cell));
	    	//System.out.println("Diameter: " + nearby_cell.getDiameter());
	    }*/

	    //System.out.println(nearby_cells.size());
/*
	for (int i=0; i<4; i++) {
	    int arg = gen.nextInt(180)+1;
	    Point vector = extractVectorFromAngle(arg);
	    if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) 
		return new Move(vector, (byte) arg);
	}*/

	// if all tries fail, just chill in place
	//return new Move(new Point(0,0), (byte)0);
/*	    
	if (memory > 0) { // follow previous direction unless it would cause a collision
	    Point vector = extractVectorFromAngle( (int)memory);
	    // check for collisions
	    if (!collides( player_cell, vector, nearby_cells, nearby_pheromes))
		return new Move(vector, memory);
	}

	// if no previous direction specified or if there was a collision, try random directions to go in until one doesn't collide
	for (int i=0; i<4; i++) {
	    int arg = gen.nextInt(180)+1;
	    Point vector = extractVectorFromAngle(arg);
	    if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) 
		return new Move(vector, (byte) arg);
	}

	// if all tries fail, just chill in place
	return new Move(new Point(0,0), (byte)0);*/
    }

    // check if moving player_cell by vector collides with any nearby cell or hostile pherome
    private boolean collides(Cell player_cell, Point vector, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
	Iterator<Cell> cell_it = nearby_cells.iterator();
	Point destination = player_cell.getPosition().move(vector);
	while (cell_it.hasNext()) {
	    Cell other = cell_it.next();
	    if ( destination.distance(other.getPosition()) < 0.5*player_cell.getDiameter() + 0.5*other.getDiameter() + 0.00011) 
		return true;
	}
	Iterator<Pherome> pherome_it = nearby_pheromes.iterator();
	while (pherome_it.hasNext()) {
	    Pherome other = pherome_it.next();
	    if (other.player != player_cell.player && destination.distance(other.getPosition()) < 0.5*player_cell.getDiameter() + 0.0001) 
		return true;
	}
	return false;
    }

    // convert an angle (in 2-deg increments) to a vector with magnitude Cell.move_dist (max allowed movement distance)
    private Point extractVectorFromAngle(int arg) {
	double theta = Math.toRadians( (double)arg );
	double dx = Cell.move_dist * Math.cos(theta);
	double dy = Cell.move_dist * Math.sin(theta);
	return new Point(dx, dy);
    }

}
