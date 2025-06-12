from lark import UnexpectedInput
import pytest
from config import ConfigParser, ExpressionTemplate, TreeToExpression, Expression, Term
from exceptions import ConfigError
import yaml
import tempfile

@pytest.fixture
def valid_config():
    return {
        "caller": {
            "name": "CustomCaller",
            "format": "GT:AD:DP:AF",
            "genotype": {
                "extract": "GT"
            },
            "vaf": {
                "extract": "AF"
            },
            "depth": {
                "extract": "DP"
            },
            "rrc": {
                "forward": {
                    "extract": "AD[0]"
                },
                "reverse": {
                    "extract": "AD[1]"
                },
                "total": {
                    "extract": "AD[0] + AD[1]"
                }
            },
            "arc": {
                "forward": {
                    "extract": "AD[0]"
                },
                "reverse": {
                    "extract": "AD[1]"
                },
                "total": {
                    "extract": "AD[0] + AD[1]"
                }
            }
        }
    }

@pytest.fixture
def config_parser():
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        return ConfigParser(f.name)

class TestConfigParser:
    @pytest.mark.schema
    def test_valid_schema(self, config_parser: ConfigParser, valid_config):
        """Test that a valid configuration passes schema validation"""
        assert config_parser.valid_schema(valid_config)
    
    @pytest.mark.schema
    def test_invalid_caller_name(self, config_parser: ConfigParser, valid_config):
        """Test that forbidden caller names are rejected"""
        invalid_config = valid_config.copy()
        invalid_config["caller"]["name"] = "BCFTools"
        assert not config_parser.valid_schema(invalid_config)
    
    @pytest.mark.schema
    def test_missing_required_fields(self, config_parser: ConfigParser, valid_config):
        """Test that missing required fields are caught"""
        invalid_config = valid_config.copy()
        del invalid_config["caller"]["format"]
        assert not config_parser.valid_schema(invalid_config)
    
    @pytest.mark.schema
    def test_invalid_format_string(self, config_parser: ConfigParser, valid_config):
        """Test that invalid format strings are caught"""
        invalid_config = valid_config.copy()
        invalid_config["caller"]["format"] = "GT;DP"
        assert not config_parser.valid_schema(invalid_config)

class TestDSLParser:
    
    @pytest.fixture
    def transformer(self):
        return TreeToExpression()
    
    @pytest.mark.dsl
    def test_simple_field(self, transformer: TreeToExpression):
        """Test parsing a simple field without metadata"""
        result = transformer.transform(ConfigParser.DSL_PARSER.parse("DP"))
        assert isinstance(result, Term)
        assert result.field == "DP"
        assert result.metadata is None
    
    @pytest.mark.dsl
    def test_field_with_index(self, transformer: TreeToExpression):
        """Test parsing a field with index metadata"""
        result = transformer.transform(ConfigParser.DSL_PARSER.parse("AD[format,0]"))
        print(result)
        assert isinstance(result, Term)
        assert result.field == "AD"
        assert result.metadata.header == "format"
        assert result.metadata.index == 0
    
    @pytest.mark.dsl
    def test_field_with_unit(self, transformer: TreeToExpression):
        """Test parsing a field with unit metadata"""
        result = transformer.transform(ConfigParser.DSL_PARSER.parse("AF[format,0,%]"))
        assert isinstance(result, Term)
        assert result.field == "AF"
        assert result.metadata.header == "format"
        assert result.metadata.index == 0
        assert result.metadata.unit == "%"
    
    @pytest.mark.dsl
    def test_simple_expression(self, transformer: TreeToExpression):
        """Test parsing a simple arithmetic expression"""
        result = transformer.transform(ConfigParser.DSL_PARSER.parse("AD[0] + AD[1]"))
        assert isinstance(result, Expression)
        assert result.operator == "+"
        assert len(result.terms) == 2
        assert all(isinstance(term, Term) for term in result.terms)
    
    @pytest.mark.dsl
    def test_complex_expression(self, transformer: TreeToExpression):
        """Test parsing a complex arithmetic expression"""
        result = transformer.transform(ConfigParser.DSL_PARSER.parse("(AD[0] + AD[1]) * 2"))
        assert isinstance(result, Expression)
        assert result.operator == "*"
        assert len(result.terms) == 2
        assert isinstance(result.terms[0], Expression)
        assert isinstance(result.terms[1], Term)
    
    @pytest.mark.dsl
    def test_invalid_expression(self, transformer: TreeToExpression):
        """Test that invalid expressions raise appropriate errors"""
        with pytest.raises(UnexpectedInput):
            transformer.transform(ConfigParser.DSL_PARSER.parse("invalid[quality]"))
    
    @pytest.mark.dsl
    def test_field_with_info_header(self, transformer: TreeToExpression):
        """Test parsing a field with info header"""
        result = transformer.transform(ConfigParser.DSL_PARSER.parse("QUAL[info]"))
        assert isinstance(result, Term)
        assert result.field == "QUAL"
        assert result.metadata.header == "info"
        assert result.metadata.index is None

class TestConfigParserIntegration:
    
    @pytest.mark.config
    def test_load_valid_config(self, config_parser, valid_config):
        """Test loading and validating a complete valid configuration"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml') as f:
            yaml.dump(valid_config, f)
            f.flush()
            config_parser.path = f.name
            result = config_parser.load()
            assert result is not None
            assert "caller" in result
    
    @pytest.mark.config
    def test_load_invalid_yaml(self, config_parser):
        """Test loading an invalid YAML file"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml') as f:
            f.write("invalid: yaml: content: [")
            f.flush()
            config_parser.path = f.name
            with pytest.raises(ConfigError):
                config_parser.load()
    
    @pytest.mark.config
    def test_load_nonexistent_file(self, config_parser):
        """Test loading a nonexistent file"""
        config_parser.path = "nonexistent.yaml"
        with pytest.raises(ConfigError):
            config_parser.load()

class TestExpressionVisitor:

    @pytest.fixture
    def transformer(self):
        return TreeToExpression()

    @pytest.mark.templating
    def test_expression_visitor(self, transformer: TreeToExpression):
        """Test the expression visitor."""
        pass


class TestExpressionTemplating:

    @pytest.mark.templating
    def test_expression_templating(self):
        """Test the expression templating system."""
        # Example format fields
        format_fields = ["DP", "AD", "AF", "GT"]
        
        # Create template handler
        template = ExpressionTemplate(format_fields)
        
        # Example 1: Simple term with index
        term_expr = {
            "type": "term",
            "field": "DP",
            "metadata": {
                "index": 0,
                "header": None,
                "unit": None
            }
        }
        
        # Example 2: Complex expression
        complex_expr = {
            "type": "expression",
            "operator": "+",
            "terms": [
                {
                    "type": "term",
                    "field": "AD",
                    "metadata": {
                        "index": 0,
                        "header": None,
                        "unit": None
                    }
                },
                {
                    "type": "term",
                    "field": "AD",
                    "metadata": {
                        "index": 1,
                        "header": None,
                        "unit": None
                    }
                }
            ]
        }
        
        # Test template formatting
        assert template.to_template(term_expr) == "DP[0]"
        assert template.to_template(complex_expr) == "(AD[0] + AD[1])"
        
        # Test Python code generation
        assert template.to_python(term_expr) == "DP[0]"
        assert template.to_python(complex_expr) == "(AD[0] + AD[1])"
            
        # Test with header and unit
        term_with_metadata = {
            "type": "term",
            "field": "AF",
            "metadata": {
                "index": 0,
                "header": "format",
                "unit": "%"
            }
        }
        
        assert template.to_template(term_with_metadata) == "AF[0,format,%]"
        assert template.to_python(term_with_metadata) == "AF[0,format,%]"