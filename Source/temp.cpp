bool Intersection(const Vertex a, const Vertex b) {
	float divide = b.position.x - a.position.x;
	if(divide < epsilon && divide > -epsilon) {
		return false;
	}
	else {
		return true;
	}
}

Vertex ComputeIntersection(const Vertex a, const Vertex b) {
	float xmax = a.w * SCREEN_WIDTH / 2;
	float divide = b.position.x - a.position.x;
	float y = a.position.y + (b.position.y - a.position.y) * (xmax - a.position.x) / divide;
	float z = a.position.z + (b.position.z - a.position.z) * (xmax - a.position.x) / divide;
	float x = xmax;

	Vertex output;
	output.position = vec3(x, y, a.position.z);
	output.w = a.w;
		
	return output;
}

void ClipSpace(vector<Vertex>& vertices) {
	for(unsigned int i = 0; i < vertices.size(); i++) {
		vertices[i].w = vertices[i].position.z / focalLength;
	}
}

bool Clip(vector<Vertex>& vertices) {

	printf("old\n");
	for(int i = 0; i < vertices.size(); i++) {
		printf("(%f,%f,%f)\n", vertices[i].position.x, vertices[i].position.y, vertices[i].position.z);
	}

	ClipSpace(vertices);

	int V = vertices.size();

	vector<Vertex> outputList = vertices;

	vector<Vertex> inputList = outputList;
	outputList.clear();

	Vertex start = inputList.back();

	for(unsigned int i = 0; i < inputList.size(); i++) {
		Vertex end = inputList[i];
		if(end.position.x < end.w * SCREEN_WIDTH / 2) {
			if(start.position.x > start.w * SCREEN_WIDTH / 2) {
				if(Intersection(start,end)) {
					outputList.push_back(ComputeIntersection(start,end));
				}
			}
			outputList.push_back(end);
		}
		else if(start.position.x < start.w * SCREEN_WIDTH / 2) {
			if(Intersection(start, end)) {
				outputList.push_back(ComputeIntersection(start, end));
			}
		}
		start.position.x = end.position.x;
		start.position.y = end.position.y;
		start.position.z = end.position.z;
		start.w = end.w;
	}

	printf("new\n");
	for(int i = 0; i < outputList.size(); i++) {
		printf("(%f,%f,%f)\n", outputList[i].position.x, outputList[i].position.y, outputList[i].position.z);
	}

	printf("outputList size: %d\n", outputList.size());
	vertices = outputList;

	return true;
}