package slather.g4;

import slather.sim.Cell;
import slather.sim.Point;
import slather.sim.GridObject;
import slather.sim.Move;
import slather.sim.Pherome;
import java.util.*;

public class Player implements slather.sim.Player {

	private Random gen;
	// weight parameters
	private final static int NUMBER_OF_RANDOM_TRY = 100;
	private final static double PHEROME_IMPORTANCE = 0.2;

	// range for tail length
	private final static int LOWEST_LEN_OF_TAIL = 10;

	// angle
	private final static int ANGLE_RANGE = 360;
	private final static int SCALE = 3; // establish mapping from ANGLE_RANGE to byte so that every arg is in range [0, 120)
	private int our_angle_range;
	
	// input variable
	private int tail;
	private double visible_distance;
	private double side;

	public void init(double d, int t, int side_length) {
		this.gen = new Random();
		this.visible_distance = d;
		this.tail = t;
		this.side = side_length;
		this.our_angle_range = Player.ANGLE_RANGE / Player.SCALE;
	}

	public Move play(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
		if (player_cell.getDiameter() >= 2) // reproduce whenever possible
			return new Move(true, (byte) -128, (byte) -128);

		// care bout your children, dude!

		// strategies choosen branch
		int arg = 0;
		int spin_sep = get_spin_sep(this.tail, this.visible_distance); 
		int detector_sep = get_detector_sep(this.tail, this.visible_distance);


		if (memory == -128 && nearby_cells.size() > 0) {
			arg = getOppositeDirection(player_cell, nearby_cells);
		} else if (nearby_cells.size() == 0) {
			arg = memory;
		} else if (isCrowded(player_cell, nearby_cells, nearby_pheromes, 2) == true) {

			int tmp_arg = 0;
			if (memory < 0 && memory >= -120) {// spin in opposite direction
				tmp_arg = spin(player_cell, (byte)(-memory - 1) , nearby_cells, nearby_pheromes, spin_sep, false);
				//arg = spin(player_cell, -memory+1, nearby_cells, nearby_pheromes, spin_sep, false);
			} else {
				tmp_arg = spin(player_cell, memory, nearby_cells, nearby_pheromes, spin_sep, true);
				//arg = spin(player_cell, memory, nearby_cells, nearby_pheromes, spin_sep, true);
			}
			
			Point tmp_vector = extractVectorFromAngle(tmp_arg);

			if (!collides(player_cell, tmp_vector, nearby_cells, nearby_pheromes)) {
				if (memory < 0 && memory >= -120) {
					return new Move(tmp_vector, (byte) (-tmp_arg - 1));
				} else {
					return new Move(tmp_vector, (byte) tmp_arg);	
				}
				
			} else {
				arg = (memory + 60) % 120;
				tmp_vector = extractVectorFromAngle(arg);
				if (!collides(player_cell, tmp_vector, nearby_cells, nearby_pheromes)) {
					
					return new Move(tmp_vector, (byte) (-arg-1));
				}

			}

		} else {
			arg = detector(player_cell, memory, nearby_cells, nearby_pheromes, detector_sep);
		}


		Point vector = extractVectorFromAngle(arg);
		if (!collides(player_cell, vector, nearby_cells, nearby_pheromes)) {
			return new Move(vector, (byte) arg);
		}

		// backup strategy: random escape
		// TODO: add scale
		for (int i = 0; i < Player.NUMBER_OF_RANDOM_TRY; i++) {
			arg = gen.nextInt(this.our_angle_range) + 1;
			vector = extractVectorFromAngle(arg);
			if (!collides(player_cell, vector, nearby_cells, nearby_pheromes))
				return new Move(vector, (byte) arg);
		}
		return new Move(new Point(0, 0), (byte) 0);	
	}

	/* get input variable determined parameters */
	private int get_spin_sep(int tail, double visible_distance) {
		return Math.max(tail, Player.LOWEST_LEN_OF_TAIL);
	}
	private int get_detector_sep(int tail, double visible_distance) {
		return 4;
	}
	private double isCrowdedThreshold(int tail, double visible_distance) {
		return 2.5;
	}

	private int getOppositeDirection(Cell player_cell, Set<Cell> nearby_cells) {
		Cell nearest = null;
		double minDis = Double.MAX_VALUE;
		for (Cell c : nearby_cells) {
			if (c.distance(player_cell) < minDis) {
				minDis = c.distance(player_cell);
				nearest = c;
			}
		}

		Point pos = nearest.getPosition();
		double vx = player_cell.getPosition().x - pos.x;
		double vy = player_cell.getPosition().y - pos.y;
		int arg = (int) extractAngleFromVector(vx, vy);
		return arg; //gen.nextInt(10) - 5;
	}

	private int getOppositeDirectionP(Cell player_cell, Set<Pherome> nearby_pheromes) {
		Pherome nearest = null;
		double minDis = Double.MAX_VALUE;
		for (Pherome c : nearby_pheromes) {
			if (c.distance(player_cell) < minDis && c.player != player_cell.player) {
				minDis = c.distance(player_cell);
				nearest = c;
			}
		}

		Point pos = nearest.getPosition();
		double vx = player_cell.getPosition().x - pos.x;
		double vy = player_cell.getPosition().y - pos.y;
		int arg = (int) extractAngleFromVector(vx, vy);
		return arg + gen.nextInt(30) - 15;
	}

	private int spin(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes, int seperation) {
		return ((this.our_angle_range) / seperation + memory) % (this.our_angle_range);
	}


	private int spin(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes, int seperation, boolean flag) {
		if (flag == true)
			return	((Player.ANGLE_RANGE / seperation) / Player.SCALE + memory) % 120 ;
		else 
			return ( (memory - (Player.ANGLE_RANGE / seperation) / Player.SCALE) % 120 + 120) %120 ;
	}

































private double angleBetweenVectors(Point p1, Point p2) {
    	double x1 = p1.x, x2 = p2.x, y1 = p1.y, y2 = p2.y;
    	double dot = x1*x2 + y1*y2;    
    	double det = x1*y2 - y1*x2;  
    	double angle = Math.atan2(det, dot);
    	if(angle < 0) angle += 2*Math.PI;
    	return angle;
    }
    private Point pathBetweenTangents(Cell player_cell, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {

    	class GridObjectAnglePair{
            GridObject gridObject;
            Point angle;
            public GridObjectAnglePair(GridObject obj, Point ang) {
                angle = ang;
                gridObject = obj;
            }
        }
        
    	
        class GridObjectAnglePairComparator implements Comparator<GridObjectAnglePair> {
            @Override
            public int compare(GridObjectAnglePair a, GridObjectAnglePair b) {
                double angle1 = Math.atan2(a.angle.y, a.angle.x);
                double angle2 = Math.atan2(b.angle.y, b.angle.x);
                if(angle1==angle2) return 0;
                return angle1<angle2?-1:1;
            }   
        }
        
        List<GridObjectAnglePair> nearby_list = new ArrayList<GridObjectAnglePair>(); 
        for(GridObject cell : nearby_cells) {
            nearby_list.add(new GridObjectAnglePair(cell, getClosestDirection(player_cell.getPosition(), cell.getPosition())));
        }
        for(GridObject pherome : nearby_pheromes) {
            if(pherome.player != player_cell.player) {
                nearby_list.add(new GridObjectAnglePair(pherome, getClosestDirection(player_cell.getPosition(), pherome.getPosition())));
            }
        }
        nearby_list.sort(new GridObjectAnglePairComparator());
        if(nearby_list.size()>1) {
        	
        	//start with biggest cell, or 0 if no cells
        	int biggest_i = 0;
        	double big_diam = 0;
        	int ind = 0;
        	for(GridObjectAnglePair g : nearby_list) {
        		if(g.gridObject instanceof Cell) {
        			if(((Cell)g.gridObject).getDiameter() > big_diam) {
        				big_diam = ((Cell)g.gridObject).getDiameter();
        				biggest_i = ind;
        			}
        		}
        		++ind;
        	}
        	
        	
            double widest = -1;
            Point widest_vector = null;
            Point prev_tangent = (Point) getTangentDirections(player_cell, nearby_list.get(biggest_i).gridObject).get(1);
            int prev_i = biggest_i;
            int sz = nearby_list.size();
            for(int i = biggest_i + 1; i < nearby_list.size() + biggest_i + 1; ++i) {
            	int k = i%sz;
            	Point p0 = nearby_list.get(prev_i).gridObject.getPosition();
            	Point p1 = nearby_list.get(k).gridObject.getPosition();
            	
            	Point current_tangent1 = (Point) getTangentDirections(player_cell, nearby_list.get(k).gridObject).get(0);
            	Point current_tangent2 = (Point) getTangentDirections(player_cell, nearby_list.get(k).gridObject).get(1);
//            	if(angleBetweenVectors(current_tangent1, current_tangent2) > Math.PI) System.out.println("No.");
//            	if(angleBetweenVectors(current_tangent1, current_tangent2) <= Math.PI) System.out.println("Yes");
                double angle = angleBetweenVectors(prev_tangent, current_tangent1);
                
                double center_angle = 
                		angleBetweenVectors(getClosestDirection(player_cell.getPosition(), p0),
                				getClosestDirection(player_cell.getPosition(), p1));
                //if overlap occurs
                if(center_angle < Math.PI && angle > Math.PI) {
                	double angle_temp = angleBetweenVectors(current_tangent2, prev_tangent);
                	if(angle_temp < Math.PI) {
                		prev_i = k;
                		prev_tangent = current_tangent2;
                	}
                //no overlap
                } else {
                	if(widest < angle) {
                		widest = angle;
                		widest_vector = prev_tangent;
                	}
                	prev_tangent = current_tangent2;
                	prev_i = k;
                }
            }
            if(widest_vector==null) {
            	//This should also return empty...
            	return pathOfLeastResistance(player_cell, nearby_cells, nearby_pheromes);
            }
            Point p2 = rotate_counter_clockwise(widest_vector, widest/2);
            if(collides(player_cell, p2, nearby_cells, nearby_pheromes)) {
            	
            	return pathOfLeastResistance(player_cell, nearby_cells, nearby_pheromes);
            }
            
            //return new Move(p3, memory);
            return p2;
        } else if(nearby_list.size() == 1) {
            return new Point(-nearby_list.get(0).angle.x,-nearby_list.get(0).angle.y);
        }
        return new Point(0,0);
    }
    /*
     * Angle in Radian
     */
    Point rotate_counter_clockwise(Point vector, double angle) {
		double newx, newy,x,y;
		x = vector.x;
		y = vector.y;
		newx = x*Math.cos(angle) - y*Math.sin(angle);
		newy = y*Math.cos(angle) + x*Math.sin(angle);
		return new Point(newx, newy);
	}
    
    
    /*
     * Returns in direction player_cell --> other
     * Second one is rotated more counter clockwise than first
     */
    private List<Point> getTangentDirections(Cell player_cell, GridObject other) {
    	List<Point> out = new ArrayList<Point>();
    	double r1 = player_cell.getDiameter()/2;
    	double r2 = 0;
    	if(other instanceof Cell) {
    		r2 = ((Cell)other).getDiameter()/2;
    	}
    	
    	double d = getDistance(player_cell.getPosition(), other.getPosition());
    	double theta2 = Math.asin((r1+r2)/d);
    	Point center_to_center = new Point(other.getPosition().x-player_cell.getPosition().x,
    			other.getPosition().y-player_cell.getPosition().y);
    	center_to_center = getUnitVector(center_to_center);
    	out.add(rotate_counter_clockwise(center_to_center, -theta2));
    	out.add(rotate_counter_clockwise(center_to_center, theta2));
    	
    	
    	return out;
    }

    private double getDistanceDirect(Point first, Point second) {
		double dist_square = (first.x - second.x)*(first.x - second.x) + (first.y - second.y)*(first.y - second.y);
    	double dist = Math.sqrt(dist_square);
    	return dist;
	}

    private double getDistance(Point first, Point second) {
    	double x = second.x;
    	double y = second.y;
    	double dist = 100;
    	for(int area_x = -1; area_x <= 1; area_x ++) {
    		for(int area_y = -1; area_y <= 1; area_y ++) {
    			x = second.x + area_x*100;
    			y = second.y + area_y*100;
    			double d = getDistanceDirect(first, new Point(x,y));
    			if( dist > d) {
    				dist = d;
    			}
    		}
    	}
    	//System.out.println("distance to " + second.x + " " + second.y +"= " + dist);
    	return dist;
    }
    //of first to second
    private Point getClosestDirection(Point first, Point second) {
    	//System.out.println(first.x +" " + first.y +" " + second.x +" "+ second.y);
    	double x = second.x;
    	double y = second.y;
    	double dist = 100;
    	Point best = null;
    	for(int area_x = -1; area_x <= 1; area_x ++) {
    		for(int area_y = -1; area_y <= 1; area_y ++) {
    			x = second.x + area_x*100;
    			y = second.y + area_y*100;
    			double d = getDistanceDirect(first, new Point(x ,y ));
    			if( dist > d) {
    				dist = d;
    				best = new Point(x-first.x,y-first.y);
    			}
    		}
    	}
    	return getUnitVector(best);
    }
    
    private Point getUnitVector(Point point) {
    	double x = point.x, y = point.y;
    	double norm = Math.hypot(x, y);
    	x /= norm;
    	y /= norm;
		
    	return new Point(x, y);
    }

	private Point pathOfLeastResistance(Cell player_cell, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {
    	class GridObjectAnglePair <GridObjectAnglePair>{
            GridObject gridObject;
            Point angle;
            public GridObjectAnglePair(GridObject obj, Point ang) {
                angle = ang;
                gridObject = obj;
            }
        }
        
        class GridObjectAnglePairComparator implements Comparator<GridObjectAnglePair> {
            @Override
            public int compare(GridObjectAnglePair a, GridObjectAnglePair b) {
                double angle1 = Math.atan2(a.angle.y, a.angle.x);
                double angle2 = Math.atan2(b.angle.y, b.angle.x);
                if(angle1==angle2) return 0;
                return angle1<angle2?-1:1;
            }   
        }
        
        List<GridObjectAnglePair> nearby_list = new ArrayList<GridObjectAnglePair>(); 
        for(GridObject cell : nearby_cells) {
        	if(getDistance(cell.getPosition(), player_cell.getPosition()) < 3)
            nearby_list.add(new GridObjectAnglePair(cell, getClosestDirection(player_cell.getPosition(), cell.getPosition())));
        }
        for(GridObject pherome : nearby_pheromes) {
            if(pherome.player != player_cell.player 
            		&& (getDistance(pherome.getPosition(), player_cell.getPosition()) < 3)) {
                nearby_list.add(new GridObjectAnglePair(pherome, getClosestDirection(player_cell.getPosition(), pherome.getPosition())));
            }
        }
        nearby_list.sort(new GridObjectAnglePairComparator());
        if(nearby_list.size()>1) {
            double widest = -1;
            int widest_index = -1;
            for(int i = 0; i < nearby_list.size(); ++i) {
            	Point p0 = nearby_list.get(i-1<0?nearby_list.size()-1:i-1).gridObject.getPosition();
            	Point p1 = nearby_list.get(i).gridObject.getPosition();
                double angle = 
                		angleBetweenVectors(getClosestDirection(player_cell.getPosition(), p0),
                				getClosestDirection(player_cell.getPosition(), p1));
                if( widest < angle ) {
                    widest = angle;
                    widest_index = i;
                }
            }
            Point p1 = nearby_list.get(widest_index).angle;
            Point p2 = nearby_list.get(widest_index-1<0?nearby_list.size()-1:widest_index-1).angle;
            
            p2 = rotate_counter_clockwise(p2, widest/2);
            return p2;
        } else if(nearby_list.size() == 1) {
            return new Point(-nearby_list.get(0).angle.x,-nearby_list.get(0).angle.y);
        }
        return new Point(0,0);
    }




































	private int findLargestGap(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes) {

		int count_pherome = nearby_pheromes.size();
		for (Pherome nearby_pherome : nearby_pheromes) {
			if (nearby_pherome.player == player_cell.player) count_pherome --;
		}
		if (count_pherome == 1 && nearby_cells.size() == 0) {
			return getOppositeDirectionP(player_cell, nearby_pheromes);
		} else if (count_pherome == 0 && nearby_cells.size() == 1) {
			return getOppositeDirection(player_cell, nearby_cells);
		}

		PriorityQueue<int[]> pq = new PriorityQueue<>((a, b) -> b[0] - a[0]);
		List<Integer> l = new ArrayList<>();

		for (Cell nearby_cell : nearby_cells) {
			double dy = nearby_cell.getPosition().y - player_cell.getPosition().y;
			double dx = nearby_cell.getPosition().x - player_cell.getPosition().x;
			l.add((int) extractAngleFromVector(dx, dy));
		}
		for (Pherome nearby_pherome : nearby_pheromes) {
			if (nearby_pherome.player == player_cell.player) continue;
			double dy = nearby_pherome.getPosition().y - player_cell.getPosition().y;
			double dx = nearby_pherome.getPosition().x - player_cell.getPosition().x;
			l.add((int) extractAngleFromVector(dx, dy));
		}

		Collections.sort(l);
		for (int i=0 ; i<l.size() ; ++i) {
			int gap = 0;
			int direction = 0;

			if (i == l.size()-1) {
				gap = l.get(0) - l.get(i) + 120;
				direction = (l.get(0) + l.get(i))/2;
			} else {
				gap = l.get(i+1) - l.get(i);
				direction = (l.get(i+1) + l.get(i))/2;
			}
			if (gap > 60) direction = (direction + 60) % 120;
			
			int[] pair = {gap, direction};
			pq.offer(pair);
		}

		int size = pq.size();

		return pq.poll()[1];
	}

	private int densityRadar(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes, int seperation) {
		double[] direction = new double[seperation];

		for (Cell nearby_cell : nearby_cells) {
			double m_x = nearby_cell.getPosition().x - player_cell.getPosition().x;
			double m_y = nearby_cell.getPosition().y - player_cell.getPosition().y;
			int index = (int) (extractRatioFromVector(m_x, m_y) * seperation);
			direction[index] += trans1(nearby_cell.distance(player_cell));
		}

		for (Pherome nearby_pherome : nearby_pheromes) {
			if (nearby_pherome.player == player_cell.player) continue;
			double m_x = nearby_pherome.getPosition().x - player_cell.getPosition().x;
			double m_y = nearby_pherome.getPosition().y - player_cell.getPosition().y;
			int index = (int) (extractRatioFromVector(m_x, m_y) * seperation);
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

		return (index * this.our_angle_range / seperation + gen.nextInt(this.our_angle_range / seperation)) % this.our_angle_range;
	}

	private int detector(Cell player_cell, byte memory, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes, int seperation) {
		if (this.visible_distance < 3) {
			// TODO: remove hardcode
			// return findLargestGap(player_cell, memory, nearby_cells, nearby_pheromes);
			Point p = pathOfLeastResistance(player_cell, nearby_cells, nearby_pheromes);
			return (int) extractAngleFromVector(p.x, p.y);
		} else {
			return densityRadar(player_cell, memory, nearby_cells, nearby_pheromes, seperation);
		}
	}


	/* use a number to reflect how dense a cell's neighborhood is */
	private double density(Cell player_cell, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes, double d_filter) {
		double weightSum = 0;
		for (Cell nearby_cell : nearby_cells) {
			if (nearby_cell.distance(player_cell) < d_filter)
				weightSum += trans1(nearby_cell.distance(player_cell));
		}

		for (Pherome nearby_pherome : nearby_pheromes) {
			if (nearby_pherome.player != player_cell.player && nearby_pherome.distance(player_cell) < d_filter) {
				weightSum += Player.PHEROME_IMPORTANCE * trans1(nearby_pherome.distance(player_cell));
			}
		}
		return weightSum;
	}
	private boolean isCrowded(Cell player_cell, Set<Cell> nearby_cells, Set<Pherome> nearby_pheromes, double d_filter) {
		return density(player_cell, nearby_cells, nearby_pheromes, d_filter) >= isCrowdedThreshold(this.tail, this.visible_distance);
	}

	
	/* check if a position is prohibited */
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

	/* Transform function to transform a number to different scale */
	public double trans1(double a) {
		return 1.0/(a+0.5);
	}
	public double trans2(double a) {
		return 1/Math.pow((0.01+a),2);
	}

	/* transform between angle and unit vector */
	// TODO: add scale
	private Point extractVectorFromAngle(int arg) {
		double theta = Math.toRadians((double) arg * Player.SCALE);
		double dx = Cell.move_dist * Math.cos(theta);
		double dy = Cell.move_dist * Math.sin(theta);
		return new Point(dx, dy);
	}
	private double findLeastAbs(double x1, double x2, double x3) {
		if (Math.abs(x1) < Math.abs(x2)) return Math.abs(x1) < Math.abs(x3) ? x1 : x3;
		else return Math.abs(x2) < Math.abs(x3) ? x2 : x3;
	}
	private double extractRatioFromVector(double dx, double dy) {
		// dx_ and dy_ are used to find "wraps around" distance 
		double dx_ = dx;//findLeastAbs(dx, dx+this.side, dx-this.side);
		double dy_ = dy;//findLeastAbs(dy, dy+this.side, dy-this.side);
		double theta = Math.atan(dy_/dx_);
		double ratio = (theta / (2 * Math.PI) + (dx_ > 0 ? (dy_ < 0 ? 1 : 0) : 0.5));
		return ratio;
	}
	private double extractAngleFromVector(double dx, double dy) {
		return extractRatioFromVector(dx, dy) * this.our_angle_range;
	}

}
