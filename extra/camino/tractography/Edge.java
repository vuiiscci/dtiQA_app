package tractography;

public final class Edge {

	private int startNode;
	private int endNode;
	private int count;
	private double length;
	private double weight;

	public Edge(int aStartNode, int aEndNode,int aCount, double aLength, double aWeight)
	{
		startNode = aStartNode;
		endNode = aEndNode;
		count = aCount;
		length = aLength;
		weight = aWeight;
	}
	
	public double getLength()
	{
		return this.length;
	}
	
	public void setLength(double aLength)
	{
		this.length = aLength;
	}
	
	public int getStartNode()
	{
		return this.startNode;		
	}
	
	public void setStartNode(int aStartNode)
	{
		this.startNode = aStartNode;		
	}
	
	public int getEndNode()
	{
		return this.endNode;
	}
	
	public void setEndNode(int aEndNode)
	{
		this.endNode = aEndNode;		
	}
	
	public int getCount()
	{
		return this.count;
	}
	
	public void setCount(int aCount)
	{
		this.count = aCount;		
	}
	
	public double getWeight()
	{
		return this.weight;
	}
	
	public void setWeight(double aWeight)
	{
		this.weight = aWeight;		
	}
}
