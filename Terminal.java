public class Terminal {
	int width;
	int height;

	Terminal(){
		// tput clear
		System.out.print("[H[J");
		width=80;
		height=20;
	}

	protected void finalize() throws Throwable {
		System.out.print("Terminal says - 'Bye!'\n");
		super.finalize();
	}
}
