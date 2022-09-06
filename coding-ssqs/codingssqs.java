import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

public class codingssqs {
	public static void main ( String args[] )
	{
		double c_i = 0;
		int      i = 0;
		double s_i, d_i;
		double a_i = 0, d_sum = 0, w_sum = 0, s_sum = 0;

		try( Scanner data = new Scanner( new File( args.length > 0 ? args[0] : "/dev/null" ))  ) {
			while( data.hasNext() ) {
				a_i = data.nextFloat();
				if( ! data.hasNext() ) {
					break;
				}
				s_i = data.nextFloat();

				++i;
				if( a_i < c_i ) {
					d_i = c_i - a_i;
					d_sum += d_i;
				} else {
					d_i = 0;
				}
				c_i = a_i + d_i + s_i;

				s_sum += s_i;
				w_sum += d_i + s_i;
			}
			data.close();

			if( i > 0 ) {
				System.out.format( "jobs                 %d\n", i );
				System.out.format( "average interarrival %f\n", (a_i/i) );
				System.out.format( "average delay        %f\n", (d_sum/i) );
				System.out.format( "average wait         %f\n", (w_sum/i) );
				System.out.format( "average service      %f\n", (s_sum/i) );
			}
			System.exit(0);
		} catch ( FileNotFoundException e ) {
			System.err.println( "Error opening data file." );
		} catch ( Exception e ) {
			System.err.println( e.getClass() );
		}
		System.exit(1);

	}
}

