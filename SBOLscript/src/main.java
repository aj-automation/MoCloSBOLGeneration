import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URI;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Set;

import javax.xml.namespace.QName;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.SuffixFileFilter;
import org.sbolstandard.core2.AccessType;
import org.sbolstandard.core2.ComponentDefinition;
import org.sbolstandard.core2.GenericTopLevel;
import org.sbolstandard.core2.OrientationType;
import org.sbolstandard.core2.SBOLConversionException;
import org.sbolstandard.core2.SBOLDocument;
import org.sbolstandard.core2.SBOLValidate;
import org.sbolstandard.core2.SBOLValidationException;
import org.sbolstandard.core2.Sequence;
import org.sbolstandard.core2.SequenceAnnotation;
import org.sbolstandard.core2.SequenceOntology;

public class main {
	
	public static class DNAPart{
		String name;
		String seq;
		String type;
		boolean revComplement;
		ArrayList<DNAPart> components = null;
		int start;
		int end;
		public DNAPart(String name, int start, int end, String sequence, String type, boolean revComplement) {
			this.name = name;
			this.seq = sequence;
			this.type = type;
			this.revComplement = revComplement;
			this.start = start;
			this.end = end;
		}
		public DNAPart(String name, int start, int end, String sequence, String type, boolean revComplement, ArrayList<DNAPart> components) {
			this.name = name;
			this.seq = sequence;
			this.type = type;
			this.revComplement = revComplement;
			this.components = components;
			this.start = start;
			this.end = end;
		}
		public void addSubComponents(ArrayList<DNAPart> components) {
			this.components = components;
		}
		public String toString() {
			return "("+this.name+","+this.start+","+this.end+","+this.type+","+this.revComplement+","+this.seq+")";
		}
	}
	
	static String version = "1";
	static String so = "http://identifiers.org/so/";
	static String provNS = "http://www.w3.org/ns/prov#";
	static String dcNS = "http://purl.org/dc/elements/1.1/";
	static String dcTermsNS = "http://purl.org/dc/terms/";
	static String celloNS = "http://cellocad.org/Terms/cello#";

	static URI activityURI;
	static String createdDate;

	
	private static URI getRole(String type) {
		if (type.equals("ribozyme")) {
	        return URI.create(so + "SO:0001977");
	    } else if (type.equals("source")){
	    	return SequenceOntology.ENGINEERED_REGION;
	    } else if (type.equals("primer_bind") || type.equals("primer")) {
	    	return URI.create(so + "SO:0005850");
	    } else if (type.equals("rep_origin")) {
	    	return URI.create(so + "SO:0000296");
	    } else if (type.equals("plasmid_vector")) {
	    	return URI.create(so + "SO:0000755");
	    } else if (type.equals("scar") || type.contains("scar")) {
	        return URI.create(so + "SO:0001953");
	    } else if (type.equals("cds") || type.equals("repressor")) {
	        return URI.create(so + "SO:0000316");
	    } else if (type.equals("promoter")) {
	        return URI.create(so + "SO:0000167");
	    } else if (type.equals("rbs")) {
	        return URI.create(so + "SO:0000139");
	    } else if (type.equals("cassette")) {
	    	return URI.create(so + "SO:0005853");
	    } else if (type.equals("terminator")) {
	        return URI.create(so + "SO:0000141");
	    } else if (type.equals("insulator")) {
	        return URI.create(so + "SO:0000627");
	    } else if (type.equals("grna")) {
	        return URI.create(so + "SO:0001264");
	    } else if (type.equals("UAS")) {
	        return URI.create(so + "SO:0001678");
	    } else if (type.equals("LinkingSequence")) {
	        return URI.create(so + "SO:0001678");
	    } else if (type.equals("TataBox")) {
	        return URI.create(so + "SO:0000174");
	    } else if (type.equals("tss")) {
	        return URI.create(so + "SO:0000315");
	    } else if (type.equals("Kozak")) {
	        return URI.create(so + "SO:0001647");
	    } else if (type.equals("Codon")) {
	        return URI.create(so + "SO:0000360");
	    } else if (type.equals("operator")) {
	        return URI.create(so + "SO:0000057");
	    } else if (type.equals("backbone")) {
	        return URI.create(so + "SO:0000755");
	    } else if (type.equals("circular")) {
	        return URI.create(so + "SO:0000988");
	    } else if (type.equals("level 2 el")){
	    	return URI.create(so + "SO:0000324");
	    } else {
	        System.err.println("Part Type " + type + " not found");
	        return null;
	    }
	}
	
	public static int isInside(int[][] arr, int[] inds) {
		int i = 0;
		for(int[] a : arr) {
			if (inds[0] >= a[0] && inds[1] <= a[1]) {
				return i;
			}
			i++;
		}
		return -1;
	}

	public static void main(String[] args) throws SBOLValidationException, IOException, SBOLConversionException {
		System.out.println("OK");
		SBOLDocument document = new SBOLDocument();
		document.setDefaultURIprefix("http://dummy.org/");
		document.setComplete(true);
		document.setCreateDefaults(true);
		
		String provNS = "http://www.w3.org/ns/prov#";
		GenericTopLevel genericTopLevel = document.createGenericTopLevel("DAMPMoClo",version, new QName(provNS , "Activity", "prov"));
		genericTopLevel.setName("DAMP MoClo Sequential Logic Library");
		Iterator<File> it = FileUtils.iterateFiles(new File("~/gbfiles/type/level0"), new SuffixFileFilter(".csv"), null);
        File f = it.next();
		while(it.hasNext()){
            System.out.println(f.getName().substring(0, f.getName().length()-4));
            DNAPart part = createLevel0Parts(f);
    		ArrayList<DNAPart> parts = part.components;
    		Iterator<DNAPart> iter = parts.iterator();
    		
    		DNAPart p = null;
    		while(iter.hasNext()) {
    			p = iter.next();
    			if(p.type == "misc_feature") {
    				if (f.getName().contains("promoter") || p.name.contains("promoter")){
    					p.type = "promoter";
    				}
    			}
    			
    		}
    		if(part!=null) {
    			addComponent(document, part);
    		}
    		f = it.next();
        }
		it = FileUtils.iterateFiles(new File("~/gbfiles/type/Level_1_NOR_Gates"), new SuffixFileFilter(".csv"), null);
        f = it.next();
		while(it.hasNext()){
            System.out.println(f.getName().substring(0, f.getName().length()-4));
            DNAPart part = createLevel0Parts(f);
    		ArrayList<DNAPart> parts = part.components;
    		Iterator<DNAPart> iter = parts.iterator();
    		
    		DNAPart p = null;
    		while(iter.hasNext()) {
    			p = iter.next();
    			if(p.type == "misc_feature") {
    				if (f.getName().contains("promoter")){
    					p.type = "promoter";
    				}
    			}
    			
    		}
    		if(part!=null) {
    			addComponent(document, part);
    		}
    		f = it.next();
        }
		it = FileUtils.iterateFiles(new File("~/gbfiles/type/Level_2_S-R_Latches"), new SuffixFileFilter(".csv"), null);
        f = it.next();
		while(it.hasNext()){
            System.out.println(f.getName().substring(0, f.getName().length()-4));
            DNAPart part = createLevel0Parts(f);
    		ArrayList<DNAPart> parts = part.components;
    		Iterator<DNAPart> iter = parts.iterator();
    		
    		DNAPart p = null;
    		while(iter.hasNext()) {
    			p = iter.next();
    			if(p.type == "misc_feature") {
    				if (f.getName().contains("promoter")){
    					p.type = "promoter";
    				}
    			}
    			
    		}
    		if(part!=null) {
    			addComponent(document, part);
    		}
    		f = it.next();
        }
		it = FileUtils.iterateFiles(new File("~/gbfiles/type/Level_1_Destination_Vector"), new SuffixFileFilter(".csv"), null);
        f = it.next();
		while(it.hasNext()){
            System.out.println(f.getName().substring(0, f.getName().length()-4));
            DNAPart part = createVector(f);
    		if(part!=null) {
    			addComponent(document, part);
    		}
    		f = it.next();
        }
		it = FileUtils.iterateFiles(new File("~/gbfiles/type/Level_2_Destination_Vectors"), new SuffixFileFilter(".csv"), null);
        f = it.next();
		while(it.hasNext()){
            System.out.println(f.getName().substring(0, f.getName().length()-4));
            DNAPart part = createVector(f);
    		if(part!=null) {
    			addComponent(document, part);
    		}
    		f = it.next();
        }
		SBOLValidate.validateSBOL(document,true,true,true);
			
		document.write("moclo1" + ".SBOL");
		// TODO Auto-generated method stub

	}
	
	public static DNAPart createVector(File f) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = br.readLine();
		String[] arr = line.split(",");
		return new DNAPart(f.getName().substring(0, f.getName().length()-4).replace('-', '_').replace('(', '_').replace(')', '_'),Integer.parseInt(arr[4]),Integer.parseInt(arr[5]),arr[1],arr[2].toLowerCase(),Boolean.parseBoolean(arr[3].toLowerCase()));
	}
	public static DNAPart createLevel0Parts(File f) throws IOException {
		ArrayList<DNAPart> l0part = new ArrayList<DNAPart>();
		
		BufferedReader br = new BufferedReader(new FileReader(f));
		String line = br.readLine();
		int[][] se = new int[1000][2];
		DNAPart wholePart = null;
		int i = 0;
		while(line != null) {
			String[] arr = line.split(",");
			if (!(arr[0].equals("source"))) {
//				se[i][0] = Integer.parseInt(arr[4]);
//				se[i][0] = Integer.parseInt(arr[5]);
				if (f.getName().contains("A1_AmtR.csv")) {
					System.out.println(arr[0]);
				}
				if (arr[0].contains("promoter") && arr[1].equals("misc_feature")) {
					arr[1] = "promoter";
				}
				l0part.add(new DNAPart(arr[0],Integer.parseInt(arr[4]),Integer.parseInt(arr[5]),arr[1],arr[2].toLowerCase(),Boolean.parseBoolean(arr[3].toLowerCase())));
				i++;
			}
			else {
				wholePart = new DNAPart(f.getName().replace('-', '_').substring(0, f.getName().length()-4),Integer.parseInt(arr[4]),Integer.parseInt(arr[5]),arr[1],arr[2].toLowerCase(),Boolean.parseBoolean(arr[3].toLowerCase()));
			}
			line = br.readLine();
		}
//		Iterator<DNAPart> iter = l0part.iterator();
//		DNAPart p = iter.next();
//		while(iter.hasNext()) {
//			Iterator<DNAPart> iter1 = l0part.iterator();
//		}
		br.close();
		wholePart.addSubComponents(l0part);
		
		
		
		return(wholePart);
		
	}
	
	public static void addComponent(SBOLDocument document, DNAPart d) throws SBOLValidationException {
		if (document.getSequence(d.name+"_sequence", version) != null || d.type.contains("primer")) {
			return;
		}
		System.out.println(d.name + " " + d.type);
		if (d.components != null) {
			Iterator<DNAPart> iter = d.components.iterator();
			DNAPart s = null;
			while (iter.hasNext()) {
				s = iter.next();
				addComponent(document,s);

			}
		}
		Sequence sequence = document.createSequence(d.name + "_sequence", version, d.seq, Sequence.IUPAC_DNA);
		sequence.setName(d.name+"_sequence");
		ComponentDefinition componentDefinition = document.createComponentDefinition(d.name, version, ComponentDefinition.DNA_REGION);
		componentDefinition.setName(d.name);
		if (d.components != null) {
			componentDefinition.addType(getRole("circular"));
		}
		//componentDefinition.addWasGeneratedBy(activityURI);
		//componentDefinition.createAnnotation(new QName(dcTermsNS,"created","dcTerms"), createdDate);

		componentDefinition.addRole(getRole(d.type));
		componentDefinition.addSequence(sequence);
		int annotationCount = 0;
		if (d.components != null) {
			Iterator<DNAPart> iter = d.components.iterator();
			DNAPart s = null;
			while (iter.hasNext()) {
				s = iter.next();

				boolean exists = (componentDefinition.getComponent(s.name) == null);
				componentDefinition.createComponent(exists ? s.name : s.name+"2", AccessType.PUBLIC, s.name, version);
				
				
				SequenceAnnotation sa = componentDefinition.createSequenceAnnotation("annotation"+annotationCount, 
						"range", s.start, s.end, s.revComplement ? OrientationType.REVERSECOMPLEMENT : OrientationType.INLINE);
				sa.setComponent(exists ? s.name : s.name+"2");
				annotationCount++;

		}
		}
		
	}

}
